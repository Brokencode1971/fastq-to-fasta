#!/usr/bin/env python3
"""
fastq_to_fasta_galaxy.py

Usage:
  python galaxy.py --galaxy-url https://usegalaxy.org --api-key XXXXXX \
      --fastq1 sample_R1.fastq.gz [--fastq2 sample_R2.fastq.gz] \
      [--history-name "DemoAssembly"] [--out assembled.fasta]

Notes:
- Requires: pip install bioblend
- Use small FASTQ files for public Galaxy (demo only).
- The script will try rnaSPAdes first, then Trinity, then other common assemblers.
"""

import os
import sys
import time
import argparse
from bioblend.galaxy import GalaxyInstance

# --- utility functions -----------------------------------------------------
def wait_for_history_ready(gi, history_id, poll_interval=5, timeout=60*60):
    """Wait until all datasets in a history are in a finished state (ok/error)."""
    start = time.time()
    while True:
        hist = gi.histories.show_history(history_id, contents=True)
        states = [d.get('state') for d in hist]
        if all(s in ('ok', 'error', 'discarded') for s in states):
            return hist
        if time.time() - start > timeout:
            raise TimeoutError("Timeout while waiting for history datasets to finish uploading/processing.")
        time.sleep(poll_interval)

def poll_history_for_dataset_with_ext(gi, history_id, exts=('fasta', 'fa', 'fna', 'fasta.gz', 'fa.gz', 'fna.gz'),
                                      name_contains=None, poll_interval=10, timeout=60*60):
    """Poll history until a dataset with one of the exts appears (and optionally name contains string)."""
    start = time.time()
    while True:
        contents = gi.histories.show_history(history_id, contents=True)
        for ds in contents:
            ext = ds.get('file_ext') or ''
            name = ds.get('name','')
            if ext in exts or any(name.lower().endswith('.' + e) for e in exts):
                if name_contains:
                    if name_contains.lower() in name.lower():
                        return ds
                else:
                    return ds
        if time.time() - start > timeout:
            raise TimeoutError("Timeout waiting for FASTA output in history.")
        time.sleep(poll_interval)

def find_tool_by_keywords(gi, keywords):
    """Search available tools for any of the keywords in tool id or name. Returns tool_id or None."""
    tools = gi.tools.get_tools()
    tool_map = {}
    for t in tools:
        tid = t.get('id','')
        name = t.get('name','')
        tool_map[tid] = name
    for kw in keywords:
        for tid, name in tool_map.items():
            if kw.lower() in tid.lower() or kw.lower() in (name or '').lower():
                return tid
    return None

def find_input_slots_for_reads(tool_info):
    """
    Inspect tool_info (as returned by gi.tools.show_tool) and try to find input parameter names
    for single-end or paired-end reads.
    Returns dict: {'paired': (left_name,right_name), 'single': single_name}
    """
    inputs = tool_info.get('inputs', [])
    single = None
    left = None
    right = None
    explicit_single = None
    explicit_fwd = None
    explicit_rev = None
    def walk_inputs(inp_list):
        nonlocal single, left, right, explicit_single, explicit_fwd, explicit_rev
        for inp in inp_list:
            name = inp.get('name') or inp.get('id') or ''
            label = inp.get('label') or inp.get('help') or ''
            type_ = inp.get('type', '')
            lname = (name or '').lower()
            llabel = (label or '').lower()
            combined = lname + ' ' + llabel
            # Capture explicit rnaspades-style keys if present (full path names when io_details=True)
            if 'single_reads' in lname:
                explicit_single = name
            if 'fwd_reads' in lname or 'forward_reads' in lname:
                explicit_fwd = name
            if 'rev_reads' in lname or 'reverse_reads' in lname:
                explicit_rev = name
            if type_ in ('data', 'data_collection', 'data_input'):
                if any(x in combined for x in ['left', 'read1', 'r1', 'forward', 'paired 1', 'paired_end_1']):
                    left = name
                elif any(x in combined for x in ['right', 'read2', 'r2', 'reverse', 'paired 2', 'paired_end_2']):
                    right = name
                elif 'reads' in combined or 'fastq' in combined or 'sequence' in combined:
                    if single is None:
                        single = name
            if 'inputs' in inp and isinstance(inp['inputs'], list):
                walk_inputs(inp['inputs'])
    walk_inputs(inputs)
    return {'paired': (left, right), 'single': single, 'explicit_single': explicit_single, 'explicit_fwd': explicit_fwd, 'explicit_rev': explicit_rev}

# --- main script -----------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Upload FASTQ to Galaxy and run a de novo assembler (best-effort).")
    parser.add_argument('--galaxy-url', required=False, help='Galaxy base URL (e.g., https://usegalaxy.org). If omitted, uses env GALAXY_URL.')
    parser.add_argument('--api-key', required=False, help='Galaxy API key. If omitted, uses env GALAXY_API_KEY.')
    parser.add_argument('--fastq1', required=True, help='Path to FASTQ file (R1 or single-end).')
    parser.add_argument('--fastq2', required=False, help='Path to FASTQ file (R2 for paired-end).')
    parser.add_argument('--history-name', default='fastq_to_fasta_history', help='History name to create on Galaxy.')
    parser.add_argument('--out', default='assembly.fasta', help='Local filename for downloaded FASTA.')
    parser.add_argument('--wait-timeout', type=int, default=60*60, help='Max seconds to wait for job completion.')
    args = parser.parse_args()

    galaxy_url = args.galaxy_url or os.environ.get('GALAXY_URL')
    api_key = args.api_key or os.environ.get('GALAXY_API_KEY')
    # REMOVED THE HARCODED API KEY - this is the critical fix
    if not galaxy_url or not api_key:
        print("ERROR: Galaxy URL and API key must be provided via --galaxy-url/--api-key or env GALAXY_URL/GALAXY_API_KEY.")
        sys.exit(1)

    if not os.path.exists(args.fastq1):
        print("ERROR: fastq1 not found:", args.fastq1); sys.exit(1)
    if args.fastq2 and not os.path.exists(args.fastq2):
        print("ERROR: fastq2 not found:", args.fastq2); sys.exit(1)

    print("Connecting to Galaxy at", galaxy_url)
    gi = GalaxyInstance(url=galaxy_url, key=api_key)

    # --- FIXED LINE: create history with named argument (was passing dict previously) ---
    print("Creating history:", args.history_name)
    history = gi.histories.create_history(name=args.history_name)
    history_id = history['id']

    # upload files
    print("Uploading FASTQ file(s)...")
    uploaded = []
    for fpath in [args.fastq1, args.fastq2] if args.fastq2 else [args.fastq1]:
        if fpath is None:
            continue
        print(" -> uploading", fpath)
        gi.tools.upload_file(fpath, history_id)
        uploaded.append(os.path.basename(fpath))

    print("Waiting for uploads to complete...")
    wait_for_history_ready(gi, history_id, poll_interval=3, timeout=args.wait_timeout//4)

    contents = gi.histories.show_history(history_id, contents=True)
    name_to_dsid = {}
    for ds in contents:
        if ds.get('deleted', False):
            continue
        ds_name = ds.get('name','')
        ds_id = ds.get('id')
        for up in uploaded:
            if up in ds_name:
                name_to_dsid[up] = ds_id

    print("Uploaded dataset IDs:", name_to_dsid)
    if len(name_to_dsid) == 0:
        print("ERROR: uploaded datasets not found in history.")
        sys.exit(1)

    candidates = ['rnaSPAdes', 'rnaspades', 'trinity', 'rnaplants', 'rna-bloom', 'rna_bloom', 'transabyss', 'oases']
    tool_id = find_tool_by_keywords(gi, candidates)
    if not tool_id:
        print("ERROR: Could not find a known assembler tool on this Galaxy instance.")
        sys.exit(1)
    print("Selected tool id:", tool_id)

    tool_info = gi.tools.show_tool(tool_id, io_details=True)
    read_slots = find_input_slots_for_reads(tool_info)
    print("Detected read input slots:", read_slots)

    tool_inputs = {}
    dsid = name_to_dsid.get(os.path.basename(args.fastq1))
    
    # Check if this is rnaSPAdes and build the correct nested structure
    if 'rnaspades' in tool_id.lower():
        print("Detected rnaSPAdes - building nested parameter structure...")
        if args.fastq2:
            # Paired-end: use "separate" mode with fwd_reads and rev_reads
            right_ds = name_to_dsid.get(os.path.basename(args.fastq2))
            if not right_ds:
                print("ERROR: Cannot find dataset id for fastq2:", args.fastq2)
                sys.exit(1)
            tool_inputs['libraries_0|files_0|file_type|type'] = 'separate'
            tool_inputs['libraries_0|files_0|file_type|fwd_reads'] = {'src': 'hda', 'id': dsid}
            tool_inputs['libraries_0|files_0|file_type|rev_reads'] = {'src': 'hda', 'id': right_ds}
        else:
            # Single-end: use "unpaired" mode with unpaired_reads
            tool_inputs['libraries_0|files_0|file_type|type'] = 'unpaired'
            tool_inputs['libraries_0|files_0|file_type|unpaired_reads'] = {'src': 'hda', 'id': dsid}
    else:
        # Fallback to original logic for other tools
        if args.fastq2:
            left_slot, right_slot = read_slots['paired']
            # Prefer explicit rnaspades keys if available
            fwd_slot = read_slots.get('explicit_fwd') or left_slot
            rev_slot = read_slots.get('explicit_rev') or right_slot
            if fwd_slot and rev_slot:
                left_ds = name_to_dsid.get(os.path.basename(args.fastq1))
                right_ds = name_to_dsid.get(os.path.basename(args.fastq2))
                if not left_ds or not right_ds:
                    print("ERROR: Cannot find dataset ids for paired files.")
                    sys.exit(1)
                tool_inputs[fwd_slot] = {'src':'hda', 'id': left_ds}
                tool_inputs[rev_slot] = {'src':'hda', 'id': right_ds}
            else:
                single_slot = read_slots['single']
                if single_slot:
                    print("Assembler expects single read input; building a dataset collection for paired upload...")
                    collection_desc = {
                        "name": "paired_collection_demo",
                        "elements": [
                            {"name": os.path.basename(args.fastq1), "src": "hda", "id": name_to_dsid[os.path.basename(args.fastq1)]},
                            {"name": os.path.basename(args.fastq2), "src": "hda", "id": name_to_dsid[os.path.basename(args.fastq2)]}
                        ]
                    }
                    coll = gi.histories.create_dataset_collection(history_id, collection_desc)
                    coll_id = coll['id']
                    tool_inputs[single_slot] = {'src':'hdca', 'id': coll_id}
                else:
                    print("ERROR: Could not find appropriate read input slots for paired-end input.")
                    sys.exit(1)
        else:
            # Prefer explicit single_reads key if present
            explicit_single = read_slots.get('explicit_single')
            single_slot = explicit_single or read_slots['single']
            if single_slot:
                tool_inputs[single_slot] = {'src':'hda', 'id': dsid}
            else:
                # Prefer a data_collection input if available (rnaSPAdes often expects a list)
                def find_collection_input(inp_list):
                    preferred = None
                    for i in inp_list:
                        type_ = i.get('type')
                        name = i.get('name') or i.get('id')
                        label = (i.get('label') or i.get('help') or '').lower()
                        if type_ == 'data_collection':
                            # Prefer names/labels mentioning reads/fastq
                            if any(k in ((name or '').lower() + ' ' + label) for k in ['read', 'fastq', 'sequence', 'rna']):
                                return name
                            if preferred is None:
                                preferred = name
                        if 'inputs' in i and isinstance(i['inputs'], list):
                            r = find_collection_input(i['inputs'])
                            if r:
                                return r
                    return preferred
                collection_slot = find_collection_input(tool_info.get('inputs', []))
                if collection_slot:
                    collection_desc = {
                        "name": "reads_list",
                        "collection_type": "list",
                        "elements": [
                            {"name": os.path.basename(args.fastq1), "src": "hda", "id": dsid}
                        ]
                    }
                    coll = gi.histories.create_dataset_collection(history_id, collection_desc)
                    coll_id = coll['id']
                    tool_inputs[collection_slot] = {'src':'hdca', 'id': coll_id}
                else:
                    # Last resort: any data input
                    def find_any_data_input(inp_list):
                        for i in inp_list:
                            if i.get('type') in ('data', 'data_input'):
                                return i.get('name') or i.get('id')
                            if 'inputs' in i:
                                r = find_any_data_input(i['inputs'])
                                if r:
                                    return r
                        return None
                    any_slot = find_any_data_input(tool_info.get('inputs', []))
                    if any_slot:
                        tool_inputs[any_slot] = {'src':'hda', 'id': dsid}
                    else:
                        print("ERROR: Could not find any suitable data input on tool.")
                        sys.exit(1)

    print("Prepared tool inputs (best-effort):", tool_inputs)
    if not tool_inputs:
        print("Assembler inputs not detected. Falling back to FASTQ→FASTA conversion (prototype mode)...")
        conv_candidates = [
            'toolshed.g2.bx.psu.edu/repos/devteam/fastq_to_fasta/fastq_to_fasta',
            'seqtk_seq',
            'seqtk_seq_v2'
        ]
        conv_tool = find_tool_by_keywords(gi, conv_candidates)
        if not conv_tool:
            print("ERROR: Could not find a FASTQ→FASTA conversion tool.")
            sys.exit(1)
        print("Selected conversion tool:", conv_tool)
        dsid = name_to_dsid.get(os.path.basename(args.fastq1))
        conv_inputs = {}
        # common input names for fastq_to_fasta or seqtk seq
        for key in ['input_fastq', 'input_file', 'fastq', 'input']:
            conv_inputs[key] = {'src': 'hda', 'id': dsid}
        # seqtk needs -a to output FASTA; many wrappers expose as 'to_fasta'
        for key in ['to_fasta', 'fasta_output', 'convert_to_fasta']:
            conv_inputs[key] = True
        print("Launching conversion tool...")
        run_resp = gi.tools.run_tool(history_id, conv_tool, conv_inputs)
        print("Tool run response:", run_resp)
        print("Waiting for conversion to finish...")
        wait_for_history_ready(gi, history_id, poll_interval=5, timeout=args.wait_timeout//2)
        print("Looking for FASTA output in history (conversion)...")
        fasta_ds = poll_history_for_dataset_with_ext(gi, history_id, exts=('fasta', 'fa', 'fna', 'fasta.gz', 'fa.gz', 'fna.gz'),
                                                    name_contains=None, poll_interval=5, timeout=900)
    else:
        print("Launching tool... this returns immediately and a job will run on the Galaxy server.")
        run_resp = gi.tools.run_tool(history_id, tool_id, tool_inputs)
        print("Tool run response:", run_resp)
        print("Waiting for tool to finish (this may take a long time depending on dataset size and instance quotas)...")
        wait_for_history_ready(gi, history_id, poll_interval=10, timeout=args.wait_timeout)
        print("Looking for FASTA output in history...")
        try:
            fasta_ds = poll_history_for_dataset_with_ext(gi, history_id, exts=('fasta', 'fa', 'fna', 'fasta.gz', 'fa.gz', 'fna.gz'),
                                                        name_contains=None, poll_interval=5, timeout=1200)
        except TimeoutError:
            contents = gi.histories.show_history(history_id, contents=True)
            print("Timed out waiting for FASTA. History contents:")
            for ds in contents:
                print(" -", ds.get('id'), ds.get('name'), ds.get('file_ext'), ds.get('state'))
            sys.exit(1)

    fasta_id = fasta_ds['id']
    fasta_name = fasta_ds.get('name')
    print("Found FASTA dataset:", fasta_id, fasta_name)

    out_path = args.out
    print("Downloading FASTA to", out_path)
    gi.datasets.download_dataset(fasta_id, out_path, use_default_filename=False)
    print("Download complete. Output saved to:", out_path)
    print("Done. History URL (view in web browser):", os.path.join(galaxy_url.rstrip('/'), 'h', history_id))

if __name__ == '__main__':
    main()
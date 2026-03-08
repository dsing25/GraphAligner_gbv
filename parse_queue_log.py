#!/usr/bin/env python3
"""
Parse FullTraceQueue.log to extract target nodes and skipFirst status.
Outputs to CSV and JSON formats.
"""

import re
import json
import csv
from collections import defaultdict

def parse_queue_log(log_file):
    """
    Parse the queue log file and extract target nodes with skipFirst status.
    Returns a list of dictionaries with node info.
    """
    # Pattern to match queue entries: [priority=X target=YYYY incoming.scoreEnd=ZZ skipFirst=0/1]
    pattern = r'\[priority=(\d+)\s+target=(\d+)\s+incoming\.scoreEnd=(\d+)\s+skipFirst=([01])\]'

    entries = []
    seen_combinations = set()  # Track unique (target, skipFirst) combinations

    with open(log_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            matches = re.findall(pattern, line)
            for match in matches:
                priority, target, score_end, skip_first = match

                entry = {
                    'line_number': line_num,
                    'target_node': int(target),
                    'skipFirst': int(skip_first),
                    'priority': int(priority),
                    'incoming_scoreEnd': int(score_end)
                }

                entries.append(entry)

                # Track unique combinations
                combo = (int(target), int(skip_first))
                seen_combinations.add(combo)

    return entries, seen_combinations

def generate_summary(entries):
    """Generate summary statistics."""
    node_stats = defaultdict(lambda: {'skipFirst_0': 0, 'skipFirst_1': 0, 'total': 0})

    for entry in entries:
        node = entry['target_node']
        skip = entry['skipFirst']

        node_stats[node]['total'] += 1
        if skip == 0:
            node_stats[node]['skipFirst_0'] += 1
        else:
            node_stats[node]['skipFirst_1'] += 1

    return dict(node_stats)

def main():
    log_file = '/data3/dsinghan/pgBench_debug/pangenomicsBench/Gbv/FullTraceQueue.log'

    print(f"Parsing {log_file}...")
    entries, unique_combos = parse_queue_log(log_file)

    print(f"\nFound {len(entries)} total queue entries")
    print(f"Found {len(unique_combos)} unique (target_node, skipFirst) combinations")

    # Generate summary
    node_stats = generate_summary(entries)
    print(f"Found {len(node_stats)} unique target nodes")

    # Save all entries to CSV
    csv_file = '/data3/dsinghan/pgBench_debug/pangenomicsBench/Gbv/queue_entries.csv'
    with open(csv_file, 'w', newline='') as f:
        if entries:
            writer = csv.DictWriter(f, fieldnames=entries[0].keys())
            writer.writeheader()
            writer.writerows(entries)
    print(f"\nSaved all entries to: {csv_file}")

    # Save unique combinations to CSV
    unique_csv = '/data3/dsinghan/pgBench_debug/pangenomicsBench/Gbv/unique_node_skipfirst.csv'
    with open(unique_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['target_node', 'skipFirst'])
        for target, skip in sorted(unique_combos):
            writer.writerow([target, skip])
    print(f"Saved unique combinations to: {unique_csv}")

    # Save summary statistics to JSON
    json_file = '/data3/dsinghan/pgBench_debug/pangenomicsBench/Gbv/node_statistics.json'
    with open(json_file, 'w') as f:
        json.dump({
            'total_entries': len(entries),
            'unique_nodes': len(node_stats),
            'unique_combinations': len(unique_combos),
            'node_stats': node_stats
        }, f, indent=2)
    print(f"Saved statistics to: {json_file}")

    # Print sample statistics
    print("\n=== Sample Node Statistics ===")
    sample_nodes = sorted(node_stats.keys())[:10]
    for node in sample_nodes:
        stats = node_stats[node]
        print(f"Node {node}: total={stats['total']}, skipFirst=0: {stats['skipFirst_0']}, skipFirst=1: {stats['skipFirst_1']}")

    if len(node_stats) > 10:
        print(f"... and {len(node_stats) - 10} more nodes")

if __name__ == '__main__':
    main()

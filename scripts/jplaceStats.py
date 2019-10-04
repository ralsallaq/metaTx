#!/usr/bin/env python

import argparse
import logging
import csv
import json
import sys
import gzip
import bz2
import os


logging.basicConfig(
    format='%(message)s',
    level=logging.INFO)
log = logging.getLogger(__name__)


# Handle py2 and py3 file differences
if sys.version_info[0] == 3:
    from io import IOBase
    file = IOBase


def open_file(path):
    if not os.path.exists(path):
        log.error("%s does not exist", path)
        return None
    # Implicit else
    suffix = path.split('.')[-1].lower().strip()
    if suffix == 'gz':
        return gzip.open(path, mode='rt')
    elif suffix == 'bz2':
        return bz2.BZ2File(path, mode='r')
    else:
        return open(path, mode='rt')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'jplace',
        metavar='XXX.jplace[.bz2|.gz]',
        help='jplace file (default is standard in)',
        nargs='?',
        default=sys.stdin,
    )
    parser.add_argument(
        '--csv',
        metavar='XXX.csv',
        help='CSV file with stats (default is standard out)',
        nargs='?',
        default=sys.stdout,
    )
    args = parser.parse_args()

    if isinstance(args.csv, file):
        csv_h = args.csv
    else:
        try:
            csv_h = open(args.csv, mode='wt')
        except Exception as e:
            log.error("Failed to open %s for writing.", args.csv)
            sys.exit(-1)
    # Implicit else we were able to open the csv file
    if isinstance(args.jplace, file):
        jplace_h = args.jplace
    else:
        jplace_h = open_file(args.jplace)

    if jplace_h is None:
        log.error("Could not open jplace file. Exiting")
        sys.exit(-1)
    # Implicit else, parse the json from the file handle
    try:
        jplace = json.load(jplace_h)
        jplace_h.close()
    except Exception as e:
        log.error("Failed to parse JSON with error {}".format(e))
        sys.exit(-1)
    # Implcit else we have something JSON can work with
    try:
        assert('fields' in jplace.keys()), "JSON Missing fields entry"
        assert('placements' in jplace.keys()), "JSON Missing placements entry"
    except AssertionError as e:
        log.error("Malformed jplace, {}".format(e))
        sys.exit(-1)
    # Implicit else we have some data...
    # Associate a field name with an index
    field_idx = {
        fn: i
        for i, fn in enumerate(jplace['fields'])
    }
    log.info("There are %d placements", len(jplace['placements']))
    writer = csv.writer(csv_h)
    writer.writerow([
        "seq_id",
        "richness",
        'min_distance'
    ])
    for placement in jplace['placements']:
        placement_sv = [s[0] for s in placement['nm']]
        placement_dl = [pl[field_idx['distal_length']] for pl in placement['p']]
        for sv in placement_sv:
            writer.writerow([
                sv,
                len(placement_dl),
                min(placement_dl),
            ])
    csv_h.close()
    sys.exit(0)


if __name__ == '__main__':
    main()

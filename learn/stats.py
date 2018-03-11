# ===-- klee-stats --------------------------------------------------------===##
# 
#                      The KLEE Symbolic Virtual Machine
# 
#  This file is distributed under the University of Illinois Open Source
#  License. See LICENSE.TXT for details.
# 
# ===----------------------------------------------------------------------===##

"""Output statistics logged by Klee."""

# use '/' to mean true division and '//' to mean floor division
from __future__ import division
from __future__ import print_function

import os
import re
import sys
import argparse

from operator import itemgetter

try:
    from tabulate import TableFormat, Line, DataRow, tabulate
except:
    print('Error: Package "tabulate" required for table formatting. '
          'Please install it using "pip" or your package manager.',
          file=sys.stderr)
    exit(1)

Legend = [
    ('Instrs', 'number of executed instructions'),
    ('Time', 'total wall time (s)'),
    ('TUser', 'total user time'),
    ('ICov', 'instruction coverage in the LLVM bitcode (%)'),
    ('BCov', 'branch coverage in the LLVM bitcode (%)'),
    ('ICount', 'total static instructions in the LLVM bitcode'),
    ('TSolver', 'time spent in the constraint solver'),
    ('States', 'number of currently active states'),
    ('Mem', 'megabytes of memory currently used'),
    ('Queries', 'number of queries issued to STP'),
    ('AvgQC', 'average number of query constructs per query'),
    ('Tcex', 'time spent in the counterexample caching code'),
    ('Tfork', 'time spent forking'),
    ('TResolve', 'time spent in object resolution'),
]

KleeTable = TableFormat(lineabove=Line("-", "-", "-", "-"),
                        linebelowheader=Line("-", "-", "-", "-"),
                        linebetweenrows=None,
                        linebelow=Line("-", "-", "-", "-"),
                        headerrow=DataRow("|", "|", "|"),
                        datarow=DataRow("|", "|", "|"),
                        padding=0,
                        with_header_hide=None)


def getLogFile(path):
    """Return the path to run.stats."""
    return os.path.join(path, 'run.stats')


class LazyEvalList:
    """Store all the lines in run.stats and eval() when needed."""

    def __init__(self, lines):
        # The first line in the records contains headers.
        self.lines = lines[1:]

    def __getitem__(self, index):
        if isinstance(self.lines[index], str):
            self.lines[index] = eval(self.lines[index])
        return self.lines[index]

    def __len__(self):
        return len(self.lines)


def getMatchedRecordIndex(records, column, target):
    """Find target from the specified column in records."""
    target = int(target)
    lo = 0
    hi = len(records) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if column(records[mid]) <= target:
            lo = mid + 1
        else:
            hi = mid
    return lo


def aggregateRecords(records):
    # index for memUsage and stateCount in run.stats
    memIndex = 6
    stateIndex = 5

    # maximum and average memory usage
    memValues = list(map(itemgetter(memIndex), records))
    maxMem = max(memValues) / 1024 / 1024
    avgMem = sum(memValues) / len(memValues) / 1024 / 1024

    # maximum and average number of states
    stateValues = list(map(itemgetter(stateIndex), records))
    maxStates = max(stateValues)
    avgStates = sum(stateValues) / len(stateValues)

    return (maxMem, avgMem, maxStates, avgStates)


def stripCommonPathPrefix(paths):
    paths = map(os.path.normpath, paths)
    paths = [p.split('/') for p in paths]
    zipped = zip(*paths)
    i = 0
    for i, elts in enumerate(zipped):
        if len(set(elts)) > 1:
            break
    return ['/'.join(p[i:]) for p in paths]


def getKeyIndex(key, labels):
    """Get the index of the specified key in labels."""

    def normalizeKey(key):
        return re.split('\W', key)[0]

    for i, title in enumerate(labels):
        if normalizeKey(title) == normalizeKey(key):
            return i
    else:
        raise ValueError('invalid key: {0}'.format(key))


def getKleeOutDirs(dirs):
    kleeOutDirs = []
    for dir in dirs:
        if os.path.exists(os.path.join(dir, 'info')):
            kleeOutDirs.append(dir)
        else:
            for root, subdirs, _ in os.walk(dir):
                for d in subdirs:
                    path = os.path.join(root, d)
                    if os.path.exists(os.path.join(path, 'info')):
                        kleeOutDirs.append(path)
    return kleeOutDirs


def getLabels(pr):
    labels = ('Path', 'Instrs', 'Time(s)', 'ICov(%)', 'BCov(%)', 'ICount',
              'TSolver(%)', 'States', 'maxStates', 'avgStates', 'Mem(MB)',
              'maxMem(MB)', 'avgMem(MB)', 'Queries', 'AvgQC', 'Tcex(%)',
              'Tfork(%)')

    return labels


def getRow(record, stats, pr):
    """Compose data for the current run into a row."""
    I, BFull, BPart, BTot, T, St, Mem, QTot, QCon, \
    _, Treal, SCov, SUnc, _, Ts, Tcex, Tf, Tr = record
    maxMem, avgMem, maxStates, avgStates = stats

    # special case for straight-line code: report 100% branch coverage
    if BTot == 0:
        BFull = BTot = 1

    Mem = Mem / 1024 / 1024
    AvgQC = int(QCon / max(1, QTot))

    row = (I, Treal, 100 * SCov / (SCov + SUnc),
           100 * (2 * BFull + BPart) / (2 * BTot), SCov + SUnc,
           100 * Ts / Treal, St, maxStates, avgStates,
           Mem, maxMem, avgMem, QTot, AvgQC,
           100 * Tcex / Treal, 100 * Tf / Treal)

    return row


def main(klee_dir):
    # function for sanitizing arguments
    def isPositiveInt(value):
        try:
            value = int(value)
        except ValueError:
            raise argparse.ArgumentTypeError(
                'integer expected: {0}'.format(value))
        if value <= 0:
            raise argparse.ArgumentTypeError(
                'positive integer expected: {0}'.format(value))
        return value


    # get print controls
    pr = 'all'

    dirs = getKleeOutDirs(klee_dir)
    if len(dirs) == 0:
        print('no klee output dir found', file=sys.stderr)
        exit(1)
    # read contents from every run.stats file into LazyEvalList
    data = [LazyEvalList(list(open(getLogFile(d)))) for d in dirs]
    if len(data) > 1:
        dirs = stripCommonPathPrefix(dirs)
    # attach the stripped path
    data = list(zip(dirs, data))

    labels = getLabels(pr)
    # labels in the same order as in the run.stats file. used by --compare-by.
    # current impl needs monotonic values, so only keep the ones making sense.
    rawLabels = ('Instrs', '', '', '', '', '', '', 'Queries',
                 '', '', 'Time', 'ICov', '', '', '', '', '', '')

    # if args.compBy:
    #     # index in the record of run.stats
    #     compIndex = getKeyIndex(args.compBy, rawLabels)
    #     if args.compAt:
    #         if args.compAt == 'last':
    #             # [records][last-record][compare-by-index]
    #             refValue = min(map(lambda r: r[1][-1][compIndex], data))
    #         else:
    #             refValue = args.compAt
    #     else:
    #         refValue = data[0][1][-1][compIndex]

    # build the main body of the table
    table = []
    totRecords = []  # accumulated records
    totStats = []  # accumulated stats
    for path, records in data:
        row = [path]
        # if args.compBy:
        #     matchIndex = getMatchedRecordIndex(
        #         records, itemgetter(compIndex), refValue)
        #     stats = aggregateRecords(LazyEvalList(records[:matchIndex + 1]))
        #     totStats.append(stats)
        #     row.extend(getRow(records[matchIndex], stats, pr))
        #     totRecords.append(records[matchIndex])
        # else:
        stats = aggregateRecords(records)
        totStats.append(stats)
        row.extend(getRow(records[-1], stats, pr))
        totRecords.append(records[-1])
        table.append(row)
    # calculate the total
    totRecords = [sum(e) for e in zip(*totRecords)]
    totStats = [sum(e) for e in zip(*totStats)]
    totalRow = ['Total ({0})'.format(len(table))]
    totalRow.extend(getRow(totRecords, totStats, pr))

    # if args.sortBy:
    #     table = sorted(table, key=itemgetter(getKeyIndex(args.sortBy, labels)),
    #                    reverse=(not args.ascending))

    if len(data) > 1:
        table.append(totalRow)
    table.insert(0, labels)

    stream = tabulate(
        table, headers='firstrow',
        tablefmt=KleeTable,
        floatfmt='.{p}f'.format(p=2),
        numalign='right', stralign='center')
    # add a line separator before the total line
    if len(data) > 1:
        stream = stream.splitlines()
        stream.insert(-2, stream[-1])
        stream = '\n'.join(stream)
    return table

if __name__ == '__main__':
    main('klee-last')

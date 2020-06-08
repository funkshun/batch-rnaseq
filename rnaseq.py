#!/usr/bin/env python

""" RNAseq Tools<br>
&nbsp;&nbsp;&nbsp;&nbsp;Online tools for processing RNA msBWTs
"""

import os, sys, time
import itertools
import argparse
import pandas as pd
import glob
import locale
import numpy as np
import matplotlib.pyplot as plot
from matplotlib.colors import LinearSegmentedColormap
import MUSCython
import MUSCython.MultiStringBWTCython as msbwt

def main():
    
    # Argument Parsing
    parser = argparse.ArgumentParser()
    
    parser.add_argument('path', 
            help='Dataset Directory, Single or Batch', 
            type=str)
    parser.add_argument('-o', '--output',
            help='Location to save output if not stdout',
            nargs=1, type=str, default='.')
    parser.add_argument('-p', '--probes',
            help='Path to file containing Variant Probes',
            nargs=1, type=str, default='./VariantProbes.csv')
    parser.add_argument('-r', '--report', 
            help='Report Style', 
            nargs=1, type=str, default='stdout',
            choices=['stdout', 'txt', 'html', 'json'])
    parser.add_argument('-s', '--separate',
            help='Generate separate report file for each dataset, ignored for stdout',
            action='store_true')
    parser.add_argument('-sdp', '--sdp_lookup',
            help='Path to file containing variant table for SDP Vector',
            nargs=1, type=str, default='./SDPpositions.csv')
    parser.add_argument('-t', '--threshold', 
            help='Minimum Voting Threshold',
            nargs=1, type=int, default=30)


    args = parser.parse_args()

    # Sanitize Path
    args.path = os.path.normpath(args.path)


    # Grab directories at target
    dcontents = sorted([x for x in glob.glob(os.path.join(args.path, '*')) if os.path.isdir(x)])

    # Determine Single or Batch
    if len(dcontents) == 0:
        if ynquery('Running Single Identification. Continue?'):
            tpaths = [args.path]
        else:
            sys.exit(0)
    else:
        if ynquery(f'Running batch identification on {len(dcontents)} samples. Continue?'):
            tpaths = dcontents
        else:
            sys.exit(0)

    reports = runQueries(tpaths, args.probes, args.sdp_lookup, args.threshold, args.report[0])

    if args.report != 'stdout':
        outputReports(reports, args.output, args.report[0], args.separate)



def runQueries(tpaths, probePath, sdpPath, threshold, report_type):
    
    # Universals
    probedf = pd.read_csv(probePath)
    sdpLookup = pd.read_csv(sdpPath)['strain']
    N = len(probedf['sdp'][0])
    
    # Iterate through found datasets and produce reports
    reports = []

    for path in tpaths:

        v, ps, t = identifySample(path, probedf, N, threshold)

        if report_type == 'stdout':
            print(createTextReport(path, sdpLookup, v, ps, t))
            continue
        elif report_type == 'txt':
            reports.append(createTextReport(path, sdpLookup, v, ps, t))
        elif report_type == 'html':
            reports.append(createHTMLReport(path, sdpLookup, v, ps, t))
        elif report_type == 'json':
            reports.append(createJSONReport(path, sdpLookup, v, ps, t))
        else:
            pass

    return reports
    
def identifySample(path, probedf, N, thresh):
    
    # Initialize Accumulators
    timing = []
    votes  = np.zeros((N,N), dtype='int32')
    ps     = 0
    
    # Load Dataset
    start = time.time()
    bwt = msbwt.loadBWT(path)
    timing.append(time.time() - start)
    
    # Voting
    start = time.time()
    for i, row in probedf.iterrows():
        
        if row['sdp'].find('H') >= 0:
            continue
        
        ref = row['probe0']
        alt = row['probe1']

        refCnt = (bwt.countOccurrencesOfSeq(senc(ref)) + 
                bwt.countOccurrencesOfSeq(senc(msbwt.reverseComplement(ref))))

        altCnt = (bwt.countOccurrencesOfSeq(senc(alt)) + 
                bwt.countOccurrencesOfSeq(senc(msbwt.reverseComplement(alt))))

        if (refCnt + altCnt) > thresh:
            
            bitVec = np.array([int(v) for v in row['sdp']], dtype='int8')
            cases = np.add.outer(bitVec, bitVec)

            if refCnt > 0:
                if altCnt > 0:
                    votes += (cases == 1)
                else:
                    votes += (cases == 0)
            else:
                votes += (cases == 2)

            ps += 1

    timing.append(time.time() - start)

    del bwt

    return votes, ps, timing

def createTextReport(path, sdps, votes, ps, timing):
    
    # Create Crosses
    indices = np.unravel_index(np.argsort(-votes,axis=None), votes.shape)
    topten  = list(zip(indices[0], indices[1]))[:10]

    report = ''
    bsep = '=' * 80 + '\n'
    ssep = '-' * 60 + '\n'

    report += bsep
    report += "Report for " + os.path.basename(path) + '\n'
    report += ssep

    report += 'Timing\n'
    report += f'BWT Loaded in {timing[0]:6.3f} seconds\n'
    report += f'Voting Concluded in {timing[1]:.3f} seconds\n'
    report += ssep

    report += 'Analysis\n'

    for i, j in topten:
        if i == j:
            cross = f'{sdps[i]}'
        else:
            cross = f'{sdps[i]} x {sdps[j]}'
        
        stats = f'{votes[i,j]}/{ps} ({(100 * (votes[i,j]/ps)):.3f}%)'

        report += cross + (' ' * (60 - (len(cross) + len(stats)))) + stats + '\n'

    report += bsep

    return report

def createHTMLReport(path, sdps, votes, ps, timing, opath='.'):
    print('HTML Reporting Not Implemented')
    return

def createJSONReport(path, sdps, votes, ps, timing, opath='.'):
    print('JSON Reporting Not Implemented')
    return

def outputReports(reports, opath, rtype, sep):
    
    if rtype == 'txt':
        if sep:
            for i, report in enumerate(reports):
                with open(os.path.join(opath, 'report' + i), 'w') as f:
                    f.write(report)
        else:
            with open(os.path.join(opath, 'batch-report'), 'w') as f:
                for report in reports:
                    f.write(report)


def senc(x):
    return x.encode('ascii', 'ignore')

def ynquery(message, default='yes'):
    
    if default is None:
        message += ' (y/n) '
    elif default == 'yes':
        message += ' (Y/n) '
    elif default == 'no':
        message += ' (y/N) '
    else:
        raise ValueError('Invalid Query Default')

    resp = input(message).lower()
    if resp in ['yes', 'y']:
        return True
    elif resp in ['no', 'n']:
        return False
    elif resp == '':
        return True if default =='yes' else False
    else:
        return None

if __name__ == '__main__':
    main()

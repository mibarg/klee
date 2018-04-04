import glob
import os
import re
import pandas as pd

timestamp_re = re.compile('(?<=\[)\d{4}.*(?=])', re.DOTALL)
min_re = re.compile('(?<=LEARN: Global minimum found: )-\d+.\d+')

def log_stats(log):
    start, end, mins = None, None, set()
    
    for r in log:
        if 'LEARN: ########## learn_weights start ##########' in r:
            m = timestamp_re.findall(r)
            start = m[0] if len(m) > 0 else None
        if 'LEARN: ########## learn_weights end ##########' in r:
            m = timestamp_re.findall(r)
            end = m[0] if len(m) > 0 else None
        if 'LEARN: Global minimum found:' in r:
            m = min_re.findall(r)
            curmin = m[0] if len(m) > 0 else 0
            mins.add(float(curmin))
    
    gmin = min(mins) if len(mins) > 0 else None
    return start, end, gmin
    

def extract_log_info(dir):

    data = pd.DataFrame(columns=['input', 'method', 'start', 'end', 'val'])

    # locate in files
    logfiles = glob.glob(os.path.join(dir, 'KleeMaze*', '*', 'log.txt'))
    for l in logfiles:
        splited = l.replace(dir, '').split(os.sep)
        with open(l, 'r') as log:
            start, end, gmin = log_stats(log)
        data = data.append(ignore_index=True, 
            other={
            'input': splited[1], 
            'method': splited[2], 
            'start': start, 
            'end': end, 
            'val': gmin})
    
    return data
            
            
DIR = '/datadrive/runtime'
data = extract_log_info(DIR)
data.to_csv(os.path.join(DIR, 'analyze.csv'))


#'[1;34m[2018/03/22 00:25:46] LEARN: ########## learn_weights start ##########'
#'[1;34m[2018/03/22 06:01:20] LEARN: Global minimum found: -34.5372'
'[1;34m[2018/03/22 06:01:21] LEARN: ########## learn_weights end ##########'
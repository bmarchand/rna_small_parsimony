import os, subprocess

os.system('rm figures/*.pdf')

list_files = []
for dataset in ['random']:
    for metric in ['num_bps','num_loops']:
        for height_method in ['min']:
            subprocess.run(['python3','average_num_bps.py',dataset,metric,height_method])

os.system('pdftk figures/*.pdf cat output figures/concatenation.pdf')
os.system('xdg-open figures/concatenation.pdf')

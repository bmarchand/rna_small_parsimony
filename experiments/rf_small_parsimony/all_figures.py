import os, subprocess


list_files = []
for dataset in ['random']:
    for metric in ['num_bps']:
        for height_method in ['min']:
            subprocess.run(['python3','average_num_bps.py',dataset,metric,height_method])

os.system('pdftk figures/*.pdf cat output figures/concatenation.pdf')
os.system('xdg-open figures/concatenation.pdf')

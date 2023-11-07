rfam_seed_lines = open('Rfam.seed',encoding='ISO-8859-1').readlines()

FAM_PREFIX = '#=GF AC' 
ID_PREFIX = '#=GF ID'
SS_PREFIX = '#=GC SS_cons'

f = open('rfam_families.txt','w')

print('\t'.join(['family','id','ss_cons']), file=f)

for line in rfam_seed_lines: 
    if line.startswith(FAM_PREFIX):
        family = line.split(' ')[-1].rstrip('\n')
    if line.startswith(ID_PREFIX):
        identifier = line.split(' ')[-1].rstrip('\n')
    if line.startswith(SS_PREFIX):
        ss_cons = line.split(' ')[-1].rstrip('\n')
        print('\t'.join([family, identifier, ss_cons]), file=f)

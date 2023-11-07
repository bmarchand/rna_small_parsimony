rfam_seed_lines = open('Rfam.seed',encoding='ISO-8859-1').readlines()

FAM_PREFIX = '#=GF AC' 
ID_PREFIX = '#=GF ID'
SS_PREFIX = '#=GC SS_cons'
INSERT_LINE_PREFIX = '#=GC RF'

f = open('rfam_families.txt','w')

print('\t'.join(['family','id','ss_cons']), file=f)

for line in rfam_seed_lines: 
    if line.startswith(FAM_PREFIX):
        family = line.split(' ')[-1].rstrip('\n')
    if line.startswith(ID_PREFIX):
        identifier = line.split(' ')[-1].rstrip('\n')
    if line.startswith(SS_PREFIX):
        ss_cons = line.split(' ')[-1].rstrip('\n')
    if line.startswith(INSERT_LINE_PREFIX):
        rf_line = line.split(' ')[-1].rstrip('\n')
        final_ss_cons = ''
        for k, c in enumerate(rf_line):
            if c!='.':
                final_ss_cons += ss_cons[k]
        print('\t'.join([family, identifier, final_ss_cons]), file=f)

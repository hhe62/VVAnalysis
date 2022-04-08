import sys,pdb

fname = sys.argv[1]
limit = float(sys.argv[2])
fin = open(fname)
iteration = 0
check = False
record = False

varl = ''
chl = ''
printvar = False
printch = False

for i,line in enumerate(fin):
    
    if 'channel:  ' in line:
        chl = line
        printch = True
        continue
    if 'variable:  ' in line:
        varl = line
        printvar = True
        continue

    if 'Iteration : ' in line:
        iteration = int(line.split('Iteration : ')[1])
        if iteration ==3:
            check = True
            continue
    
    if 'Chi^2 of change ' in line and check:
        chi2 = float(line.split('Chi^2 of change ')[1])
        check = False
        if chi2 > limit:
            print('Issue in line %s'%(i+1))
            print(line)
            record = True
            continue

    if record:
        if printvar:
            print(varl)
            printvar = False
        if printch:
            print(chl)
            printch = False
        if 'Position Indicator:' in line:
            print(line)
            record = False
            

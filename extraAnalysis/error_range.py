import numpy as np
import pdb,math
import matplotlib.pyplot as plt
import ROOT as r

years = ['2016','2017','2018']
#Choose one of the three var list:

#varstr="nJets mjj dEtajj jetPt[0] jetPt[1] absjetEta[0] absjetEta[1] MassAllj Mass0j Mass1j Mass2j Mass3j Mass4j"
varstr="nJets mjj dEtajj jetPt[0] jetPt[1] absjetEta[0] absjetEta[1]"
#varstr="MassAllj Mass0j Mass1j Mass2j Mass34j"
#varstr="MassFull Mass0jFull Mass1jFull Mass2jFull Mass34jFull"

vars = varstr.split(" ")
firstVar = varstr.split(" ")[0]

testind = 13
#vars = ["nJets"]
edict = {"stat":"statistical","pu":"PU"}
totaldict = {}

errRangeDict = {}
errOrder = []

def plotVar(errup, errdn, errtypes,var,year):
    errupa = np.array(errup)
    errdna = np.array(errdn)

    ind = errupa.argsort()
    errtypes,errup,errdn = [np.take(x, ind) for x in [errtypes, errup, errdn]]

    y = np.arange(errup.size)

    fig, axes = plt.subplots(ncols=2, sharey=True)

    #for i,j in enumerate(errdn):
    #    axes[1].text(j, i, str(j), color='black', fontweight='bold')
    #for i,j in enumerate(errup):
    #    axes[0].text(0+0.05, i, str(j), color='black', fontweight='bold')

    axes[0].barh(y, errup, align='center', color='red', zorder=10)
    axes[0].set(title='Upper Error')
    axes[1].barh(y, errdn, align='center', color='blue', zorder=10)
    axes[1].set(title='Lower Error')


    axes[0].invert_xaxis()

    axes[1].set(yticks=y, yticklabels=errtypes)
    axes[0].set_xlabel(var)
    axes[1].set_xlabel("Error Portion")
    #axes[1].yaxis.tick_right()

    for ax in axes.flat:
        ax.margins(0.03)
        ax.grid(True)

    fig.tight_layout()
    #fig.subplots_adjust(wspace=0.09)
    #plt.show()
    plt.savefig(year+"FullLog/plots/"+var+"_errorplot_%s.png"%year)
    plt.close(fig)

f16 = r.TFile('2016Full.root')
f17 = r.TFile('2017Full.root')
f18 = r.TFile('2018Full.root')
ftot = r.TFile('allFull.root')
errtypesTot=[]
errupTot=[]
errdnTot=[]

for var in vars:
    areas = {}
    h16 = f16.Get("tot_"+var+"_unf")
    h17 = f17.Get("tot_"+var+"_unf")
    h18 = f18.Get("tot_"+var+"_unf")
    htot = ftot.Get("tot_"+var+"_unf")
    sysUp =ftot.Get("tot_"+var+"_totUncUp")
    sysDn =ftot.Get("tot_"+var+"_totUncDown")
    
    a16 = h16.Integral(1,h16.GetNbinsX())
    a17 = h17.Integral(1,h17.GetNbinsX())
    a18 = h18.Integral(1,h18.GetNbinsX())
    atot = htot.Integral(1,htot.GetNbinsX())
    totStat = 0.
    totUp = 0.
    totDn = 0.
    #pdb.set_trace()
    for i in range(1,htot.GetNbinsX()+1):
        totStat += abs(htot.GetBinError(i))/atot
        #print(var,'totUp Bin',i,math.sqrt(htot.GetBinError(i)**2+sysUp.GetBinContent(i)**2)/atot)
        totUp +=math.sqrt(htot.GetBinError(i)**2+sysUp.GetBinContent(i)**2)/atot
        totDn +=math.sqrt(htot.GetBinError(i)**2+sysDn.GetBinContent(i)**2)/atot

    asum = a16+a17+a18
    areas['2016']=a16
    areas['2017']=a17
    areas['2018']=a18

    for year in years:
        start = False
        errtypes = []
        errup = []
        errdn = []

        fin = open(year+"FullLog/ErrorInfo_%s.log"%var)
        for line in fin:
            if "Error Summary" in line:
                start = True
                continue
            if start:
               if "Sum portion up and down" in line:
                   continue

               etype = line.split(":")[0] 
               if etype in edict.keys():
                   etype = edict[etype]
               errtypes.append(etype)
               eup = float(line.split("PortionUp ")[1].split(" PortionDn ")[0])
               edn = float(line.split(" PortionDn ")[1])
               errup.append(eup)
               errdn.append(edn)
        
        if year == '2016':
            errtypesTot=errtypes
            errupTot = [a16*et for et in errup]
            errdnTot = [a16*et for et in errdn]

        elif year == '2017' or year == '2018':
            for i in range(0,len(errupTot)):
                #directly set total stat error and total error here, multiplied by area since they are normalized again below
                if errtypesTot[i] == "statistical":
                    errupTot[i] = totStat*asum
                    errdnTot[i] = totStat*asum
                elif errtypesTot[i] == "total":
                    errupTot[i] = totUp*asum
                    errdnTot[i] = totDn*asum
                else:
                    errupTot[i] = errupTot[i] + areas[year]*errup[i]
                    errdnTot[i] = errdnTot[i] + areas[year]*errdn[i]

        #print('Var Year type area eup edn')
        #print(var, year, errtypes[testind],areas[year],errup[testind],errdn[testind])

        #plotVar(errup, errdn, errtypes,var,year)

    
    errupTot = [et/asum for et in errupTot]
    errdnTot = [et/asum for et in errdnTot]
    #print(var, 'Total', errtypesTot[testind],asum,atot,errupTot[testind],errdnTot[testind])
    #print(errtypesTot)
    if var == firstVar:
        for key in errtypesTot:
            if not (key == 'statistical' or key =='total'):
                errRangeDict[key] = []
    
    for order,key in enumerate(errtypesTot):
        if not (key == 'statistical' or key =='total'):
            errRangeDict[key].append(errupTot[order])
            errRangeDict[key].append(errdnTot[order])

    errOrder.append(var) #one variable has two places (up/dn) in the errRangeDict list for each err type.
    #plotVar(errupTot, errdnTot, errtypesTot,var,'Total')

#print(errRangeDict.keys())
keysReorder = [ 'eEff','mEff', 'jer', 'jes', 'fake','PU','lumi','generator', 'ggZZxsec','QCD_scales','pdf','alpha_s', ]
#pdb.set_trace()
massprint = False
for key in keysReorder:
    if not massprint:
        #print("Type: %s %s - %s"%(key,100*round(min(errRangeDict[key]),4),100*round(max(errRangeDict[key]),4) ))
        print("& %s - %s \\"%(100*round(min(errRangeDict[key]),4),100*round(max(errRangeDict[key]),4) )+"%")
        #print("Type: %s & %s \\"%(key,100*round(max(errRangeDict[key]),4) )+"%")

    else:
        pstr = ""
        for i,v in enumerate(vars):
            pstr+="& %s "%(100*round(max(errRangeDict[key][2*i:2*i+1+1]),4) )
            if round(max(errRangeDict[key][2*i:2*i+1+1]),4)>0.:
                pstr+="\\%"
            if i<len(vars)-1:
                pstr+=" "
        print(pstr)
    #print("& %s \\"%(100*round(max(errRangeDict[key]),4) )+"%")
    #print("Ind of max:%s"%(errRangeDict[key].index(max(errRangeDict[key]))))

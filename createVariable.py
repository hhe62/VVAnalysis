import json

outputname = "varsFile.json"

dict = {}
#Mass variables===================================================================
varlist0 = ["Mass"]
varlist = []
varlistFull = []

for var in varlist0:
    for j in range(0,5):
        varlist.append(var+"%sj"%j)

varlist.append("MassAllj") #for comparison check with original Mass variable
varlist.append("Mass34j") #>=3j

for var in varlist0:
    for j in range(0,5):
        varlistFull.append(var+"%sj"%j+"Full")

varlistFull.append("MassFull") #for comparison check with original Mass variable
varlistFull.append("Mass34jFull") #>=3j

print(varlist)
print(varlistFull)

for var in varlist:
    dict[var] = {}
    nj = var.replace("j","").replace("Mass","")
    if var == "MassAllj":
        nj = "All"
    if var == "Mass4j":
        nj = "\\geq 4"
    if var == "Mass34j":
        nj = "\\geq 3"
    dict[var]["units"] = '[GeV]'
    if nj in ['All','0','1']:
        dict[var]["_binning"] = [100.] + [200.+50.*i for i in range(5)] + [500.,600.,800.,1000.]
    else:
        dict[var]["_binning"] = [100.,200.,400.,600.,1000.] 

    if nj in ['0','1']:    
        dict[var]["prettyVars"] = 'm_{4\\ell}' + "(%s jet)"%nj
    elif nj == "All":
        dict[var]["prettyVars"] = 'm_{4\\ell}' + "(%sJets)"%nj
    else:
        dict[var]["prettyVars"] = 'm_{4\\ell}' + "(%s jets)"%nj
    dict[var]["responseClassNames"] = 'testJet'

for var in varlistFull:
    dict[var] = {}
    nj = var.replace("jFull","").replace("Mass","")
    if var == "MassFull":
        nj = "All"
    if var == "Mass4jFull":
        nj = "\\geq 4"
    if var == "Mass34jFull":
        nj = "\\geq 3"
    dict[var]["units"] = '[GeV]'
    if nj in ['All','0','1']:
        dict[var]["_binning"] = [80.,100.,120.,130.,180.,230.,300.,450.,600.,800.,1300.]
    else:
        dict[var]["_binning"] = [80.,100.,120.,130.,180.,230.,300.,450.,600.,800.,1300.]

    if nj in ['0','1']:
        dict[var]["prettyVars"] = 'm_{4\\ell}' + "(%s jet)"%nj
    elif nj == "All":
        dict[var]["prettyVars"] = 'm_{4\\ell}' + "(%sJets)"%nj
    else:
        dict[var]["prettyVars"] = 'm_{4\\ell}' + "(%s jets)"%nj

    dict[var]["responseClassNames"] = 'testJet'
#====================================================================================    

#jet variables=======================================================================
var2="nJets"
dict[var2] = {}
dict[var2]["units"] = ''
dict[var2]["_binning"] = [0.0,1.0,2.0,3.0,4.0]
dict[var2]["prettyVars"] = 'N_{jets}'
dict[var2]["responseClassNames"] = 'testJet'

var3="mjj"
dict[var3] = {}
dict[var3]["units"] = '[GeV]'
dict[var3]["_binning"] = [0.,200.,400.,600.,1000.]#[100.+40.*i for i in range(31)]
dict[var3]["prettyVars"] = 'Dijet Mass'
dict[var3]["responseClassNames"] = 'testJet'

var4="dEtajj"
dict[var4] = {}
dict[var4]["units"] = ''
dict[var4]["_binning"] = [0.,1.2,2.4,3.6,4.7]
dict[var4]["prettyVars"] = '|#Delta#eta(j_{1}, j_{2})|'
dict[var4]["responseClassNames"] = 'testJet'

for i,var in enumerate(["jetPt[0]","jetPt[1]","absjetEta[0]","absjetEta[1]"]):
    dict[var] = {}
    if "Pt" in var:
        dict[var]["units"] = '[GeV]'
    else:
        dict[var]["units"] = ''
    if i==0:
        dict[var]["_binning"] = [30.,50.,100.,200.,300.,500.]
        dict[var]["prettyVars"] ='p_{T}^{j1}'
    if i==1:
        dict[var]["_binning"] = [30.,50.,100.,170.,300.]
        dict[var]["prettyVars"] ='p_{T}^{j2}'
    if i==2:
        dict[var]["_binning"] = [0.,1.5,2.4,3.2,4.7]
        dict[var]["prettyVars"] = '|\\eta_{j1}|'
    if i==3:
        dict[var]["_binning"] = [0.,1.5,3.,4.7]
        dict[var]["prettyVars"] = '|\\eta_{j2}|'
    
    dict[var]["responseClassNames"] = 'testJet'    

#change axis range and MC symbol location
for key in dict.keys():
    dict[key]['top_xy']=(0.2,0.87)
    dict[key]['top_size']=0.15
    dict[key]['bottom_xy']=(0.2,0.91) #top and bottom different pad settings
    dict[key]['xyP3']=(0.2,0.91)
    dict[key]['xyP4']=(0.2,0.91)
    dict[key]['bottom_size']=0.15
    dict[key]['size_P3']=0.15
    dict[key]['size_P4']=0.15
    dict[key]['ymax_fac']=1.
    dict[key]['ymin_fac']=1.
    dict[key]['ymin_fac_extra']=1.
    dict[key]['ratio_max'] =1.8
    dict[key]['ratio_min'] = 0.4
    dict[key]['ytilt_fac'] = 1.3

    if not "Mass" in key and not "Full" in key:
        bmg = 0.45
        dict[key]['top_xy']=(0.2,0.82)
        dict[key]['top_size']=0.2
        dict[key]['bottom_xy']=(0.2,0.82) #top and bottom different pad settings
        dict[key]['xyP3']=(0.2,bmg+(1-bmg)*0.82)
        dict[key]['xyP4']=(0.2,0.82)
        dict[key]['bottom_size']=0.2
        dict[key]['size_P3']=0.2*(1-bmg)
        dict[key]['size_P4']=0.15
        dict[key]['ymax_fac']=1.
        dict[key]['ymin_fac']=1.
        dict[key]['ymin_fac_extra']=1.
        dict[key]['ratio_max'] =1.8
        dict[key]['ratio_min'] = 0.4
    if "Mass" in key and not "Full" in key:
        bmg = 0.4
        dict[key]['top_xy']=(0.2,0.2)
        dict[key]['top_size']=0.2
        dict[key]['bottom_xy']=(0.2,0.2) #top and bottom different pad settings
        dict[key]['xyP3']=(0.2,0.2)
        dict[key]['xyP4']=(0.2,bmg+(1-bmg)*0.2)
        dict[key]['bottom_size']=0.2
        dict[key]['size_P3']=0.2
        dict[key]['size_P4']=0.2*(1-bmg)
        dict[key]['ratio_max'] = 1.2
        dict[key]['ratio_min'] = 0.3
        if not "All" in key:
            dict[key]['ratio_max'] = 1.8  
            dict[key]['top_xy']=(0.2,0.82)
            dict[key]['top_size']=0.2
            dict[key]['bottom_xy']=(0.2,0.82) #top and bottom different pad settings
            dict[key]['xyP3']=(0.2,0.82)
            dict[key]['xyP4']=(0.2,bmg+(1-bmg)*0.82)
            dict[key]['bottom_size']=0.2
            dict[key]['size_P3']=0.2
            dict[key]['size_P4']=0.2*(1-bmg)  
            dict[key]['ratio_min'] = 0.3
            if "0" in key:
                dict[key]['ratio_max'] = 1.6
                dict[key]['ratio_min'] = 0.2
                

#Adjust settings for individual variables
dict['jetPt[0]']['ratio_max'] = 2.99
dict['nJets']['ratio_max'] = 2.5
dict['jetPt[0]']['ratio_min'] = 0.4
dict['jetPt[0]']['ytilt_fac'] = 1.2
dict['jetPt[1]']['ytilt_fac'] = 1.2

#dict['Mass0j']['ratio_min'] = 0.3
#ict['MassAllj']['ratio_min'] = 0.3


dict['jetPt[1]']['ymax_fac'] = 1.
dict['absjetEta[1]']['ymax_fac'] = 1.3
dict['absjetEta[0]']['ymax_fac'] = 1.3
dict['dEtajj']['ymax_fac'] = 1.2
dict['jetPt[1]']['ratio_max'] = 2.99
dict['jetPt[1]']['ratio_min'] = 0.4
#dict['jetPt[1]']['top_xy'] = (0.6,0.87)
#dict['jetPt[1]']['bottom_xy'] = (0.6,0.91)
dict['jetPt[1]']['ymin_fac'] =1.3
dict['MassFull']['ymin_fac'] =1.5
dict['jetPt[1]']['ymin_fac_extra'] =0.3

#dict['Mass4j']['top_xy'] = (0.6,0.87)
#dict['Mass4j']['bottom_xy'] = (0.6,0.91)

dict['Mass1jFull']['top_xy'] = (0.4,0.87)
dict['Mass1jFull']['bottom_xy'] = (0.4,0.91)
dict['Mass2jFull']['top_xy'] = (0.37,0.87)
dict['Mass2jFull']['bottom_xy'] = (0.37,0.91)
dict['Mass34jFull']['top_xy'] = (0.37,0.87)
dict['Mass34jFull']['bottom_xy'] = (0.37,0.91)
dict['Mass4jFull']['top_xy'] = (0.5,0.87)
dict['Mass4jFull']['bottom_xy'] = (0.5,0.91)

with open(outputname,'w') as output_file:
  json.dump(dict,output_file,indent=4)

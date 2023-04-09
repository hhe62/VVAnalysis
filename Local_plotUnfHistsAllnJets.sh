#variables="pt mass zpt leppt dphiz1z2 drz1z2"
#variables="mass"
variables="nJets mjj dEtajj jetPt[0] jetPt[1] absjetEta[0] absjetEta[1] MassAllj Mass0j Mass1j Mass2j Mass34j"
#variables="jetPt[0]" #"nJets mjj dEtajj jetPt[0] jetPt[1] absjetEta[0] absjetEta[1]" #Mass0j Mass1j Mass2j Mass34j"
#variables="MassAllj Mass0j Mass1j Mass2j Mass34j"

for var in $variables;do
  echo $var
  
  #./Utilities/scripts/plotUnfolded.py -a ZZ4l2016 -s TightLeptonsWGen -l 35.9 -f ZZ4l2016 -vr ${var} --test --makeTotals --unfoldDir /afs/cern.ch/user/u/uhussain/www/ZZFullRun2/PlottingResults/ZZ4l2016/ZZSelectionsTightLeps/ANPlots/ZZ4l2016/FinalDiffDist_16Apr2020/

  #./Utilities/scripts/plotUnfolded.py -a ZZ4l2017 -s TightLeptonsWGen -l 41.5 -f ZZ4l2017 -vr ${var} --test --makeTotals --unfoldDir /afs/cern.ch/user/u/uhussain/www/ZZFullRun2/PlottingResults/ZZ4l2017/ZZSelectionsTightLeps/ANPlots/ZZ4l2017/FinalDiffDist_16Apr2020/  
  
  #./Utilities/scripts/plotUnfolded.py -a ZZ4l2018 -s TightLeptonsWGen -l 59.7 -f ZZ4l2018 -vr ${var} --test --makeTotals --unfoldDir /afs/cern.ch/user/u/uhussain/www/ZZFullRun2/PlottingResults/ZZ4l2018/ZZSelectionsTightLeps/ANPlots/ZZ4l2018/FinalDiffDist_16Apr2020/  
  
  ./Utilities/scripts/plotUnfolded.py -a ZZ4l2018 -s TightLeptonsWGen -l 137.58 -f ZZ4l2018 -vr ${var} --test --makeTotals --scaleymin 0.3 --scaleymax 1.2 --unfoldDir unfoldPlotsTest/20230409_MiNNLO_temp_allUpdatedPlots_adjusted

done

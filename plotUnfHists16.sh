#variables="pt mass zpt leppt dphiz1z2 drz1z2"
variables="nJets mjj dEtajj jetPt[0] jetPt[1] absjetEta[0] absjetEta[1] MassAllj Mass0j Mass1j Mass2j Mass3j Mass4j"
#variables="dphiz1z2"
for var in $variables;do
  echo $var
  ./Utilities/scripts/plotUnfolded.py -a ZZ4l2016 -s TightLeptonsWGen -l 35.9 -f ZZ4l2016 -vr ${var} --test --makeTotals --scaleymin 0.3 --scaleymax 1.2 --unfoldDir /afs/cern.ch/user/h/hehe/www/FullvarList_20May2021/UnfoldZZ4l2016_oldMCreg_Fullsyst20211222/
  #--plotResponse
done

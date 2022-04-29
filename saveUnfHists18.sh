#variables="pt mass zpt leppt dphiz1z2 drz1z2"
#variables="leppt"
variables="nJets mjj dEtajj jetPt[0] jetPt[1] absjetEta[0] absjetEta[1] MassAllj Mass0j Mass1j Mass2j Mass3j Mass34j Mass4j MassFull Mass0jFull Mass1jFull Mass2jFull Mass3jFull Mass34jFull Mass4jFull"
for var in $variables;do
  echo "Now processing:$var"
  ./Utilities/scripts/saveUnfolded.py -a ZZ4l2018 -s LooseLeptons -l 59.7 -f ZZ4l2018 -sf data/scaleFactorsZZ4l2018.root -ls 2018fWUnc_full -vr ${var} --test --makeTotals --plotResponse --plotDir /afs/cern.ch/user/h/hehe/www/FullvarList_20May2021/reg_matrixPlot_2018_newMCqqZZasNom_SecondRound #--diagnostic #--noNorm
  #--makeTotals --plotResponse
done

echo "Job done!================"

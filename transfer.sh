#
# File copies simulation code over to the artemis cluster. 
#
rsync -xvr *.py *.pbs gutsim stats diet-nutrient-breakdown-SSB.xlsx art:/project/RDS-FSC-microanal-RW/gutsim/.
rsync -xvr results/20181204-trf art:/project/RDS-FSC-microanal-RW/gutsim/results/.

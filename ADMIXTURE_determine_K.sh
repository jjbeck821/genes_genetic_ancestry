#########################
# ADMIXTURE determine K #
#########################

# INPUT DATA VARIABLES
GTBOTHREF=1KGnGONL_PC_SNPs_REF

# Initial prep
plink2 --bfile $GTBOTHREF --maf 0.01 --indep-pairwise 250kb 1 0.50 --out Re-PruneBothRef
plink2 --bfile $GTBOTHREF --extract Re-PruneBothRef.prune.in --make-bed --out K-DET_1KGONL
mv K-DET_1KGONL.fam K-DET_1KGONL.fam_org
awk '{print $1" "$2" 0 0 2 -9"}' K-DET_1KGONL.fam_org | tr ' ' '\t' > K-DET_1KGONL.fam

# Run cross-validation procedure in ADMIXTURE by varying K from 3 to 27
for K in {3..27}; \
    do ./admixture -j4 K-DET_1KGONL.bed --cv $K | tee log${K}.out; done

import os
import pandas as pd
from rpy2.robjects import r, pandas2ri
import statsmodels.api as sm
# pandas2ri.activate()
os.environ["PATH"] = "/share/apps/R-3.4.2/bin:" + os.environ["PATH"]
os.environ["LD_LIBRARY_PATH"] = "/share/apps/R-3.4.2/lib64/R/lib:" + os.environ["LD_LIBRARY_PATH"]


r['load']("22.geno.RData")

batch = pandas2ri.ri2py_dataframe(r['batch'])
geno = pandas2ri.ri2py_dataframe(r("as.data.frame(geno)"))
geno = geno.astype(int)
logSigMat = pandas2ri.ri2py_dataframe(r['logSigMat'])
info = pandas2ri.ri2py_dataframe(r['info'])

data = pd.concat([logSigMat, batch, geno], axis=1)

sig = logSigMat.columns.tolist()[0]
snp = geno.columns.tolist()[0]

name1 = [sig, snp]
name1.extend(batch.columns.tolist()[14:])
name2 = [sig, snp]
name2.extend(batch.columns.tolist()[:14])

d = data.loc[:, name1]

d = d.rename(columns={d.columns[0]: "y", d.columns[1]: "x"})
d.insert(0, "group", d.TSS.str.cat(d.plate, sep=","))

model = sm.MixedLM.from_formula("y ~ x + C(cancer) - 1", data=d.dropna(), groups='group')

result = model.fit()

print(result.summary())

aic = 2 * (len(result.params) + 1) - 2 * result.llf

print(aic)


r['library']("nlme")

for sig in logSigMat.columns.tolist():
    for()

sample_list = [list(x) for x in r['sample_list']]
snp_list = [list(x) for x in r['snp_list']]


r['geno_list'].rx2(2)


# robjects.r['load']("/home/jintao/Projects/TCGA/hg19/pmSignatre/SNV10/plus_minus.RData")
# sigMat = robjects.r["pmSigMat_plus"]
# sigMat = sigMat[sigMat > 1e-11]
# sigMat = sigMat[sigMat.BG > 0]
# logSigMat = np.log(sigMat.iloc[:, :-1].div(sigMat.BG, axis=0))
# sampleInt = logSigMat.index.intersection(robjects.r['mem'].index)

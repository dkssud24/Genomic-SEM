library(GenomicSEM)
library(data.table)
library(dplyr)
files <- c("GCST90436048.tsv.input.GenomicSEM","GCST90475909.tsv.input.GenomicSEM","GCST90436053.tsv.input.GenomicSEM")
trait.names <- c("Dizziness","Hearlingloss", "Tinnitus" )
N <- c(18235, 399422, 2057)
sample.prev <- c(0.01, 0.56, 0.001271)
population.prev <- sample.prev
# 2. munging (LDSC용 summary statistics 생성)
hm3 <- "/BiO/hae/000004_ldsc/ldsc/w_hm3.snplist"
munge(files = files,
      hm3 = hm3,
      trait.names = trait.names,
      N = N)
# 3. LDSC 실행
traits <- paste0(trait.names, ".sumstats.gz")
ld <- "/BiO/hae/000004_ldsc/ldsc/eur_w_ld_chr/"
wld <- ld
ldsc_result <- ldsc(
  traits = traits,
  sample.prev = sample.prev,
  population.prev = population.prev,
  ld = ld,
  wld = wld,
  trait.names = trait.names
)
save(ldsc_result, file = "LDSC_250609.RData")
model <- "F1 =~   Dizziness+ Hearlingloss + Tinnitus"
result <- usermodel(covstruc = ldsc_result, model = model, std.lv = TRUE)
print(result$results)
print(result$modelfit)
ref <- '/BiO/hae/000006_ref_1000G/ref.freq.frq'
se.logit <- rep(FALSE, 3)
linprob <- rep(FALSE, 3)
SNPs <- sumstats(
  files = files,
  ref = ref ,
  trait.names = trait.names,
  se.logit = se.logit,
  linprob = linprob,
  N = N
)

# 6. Genomic SEM GWAS (userGWAS)
model <- "
F1 =~ Hearlingloss + Tinnitus+  Dizziness
F1 ~ SNP
"
GWAS_result <- userGWAS(
  covstruc = ldsc_result,
  SNPs = SNPs,
  model = model,
  sub = "F1~SNP",
  parallel = TRUE,
  cores = 20 
)


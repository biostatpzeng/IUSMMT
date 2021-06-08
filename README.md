# IUSMMT: survival mediation analysis of gene expression with multiple DNA methylation exposures and its application to cancers of TCGA
# Introduction
Intersection-union survival mixture-adjusted mediation test (IUSMMT) is a R procedure for examining whether a set of gene-based methylation loci affects cancer survival through gene expression under the framework of mixed models. Mediation effects are determined by two separate tests: one for the association between methylations and the expression, the other for the association between the expression and the survival outcome conditional on methylations. IUSMMT effectively combines the evidence of the two tests and infers the emergence of mediation effect by fitting an empirical three-component mixture null distribution.

Specifically, let T be a n by 2 vector of survival outcome on n individuals (including the survival time t and the survival status d), X is a n by p matrix for covariates, M is a n by K matrix for genotypes of SNPs for a genetic region (i.e. gene), G is a n by 1 vector of gene expression. For a gene under investigation, the associations between DNA methylations, expression, or clinical covariates (X) and the survival outcome can be determined within the framework of mediation analysis through the following procedures:

#### Step 1: Cox linear mixed-effects model testing for the total effect of methylations on the survival outcome ####

<a href="https://www.codecogs.com/eqnedit.php?latex=\gamma^{DE}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\gamma^{DE}" title="\gamma^{DE}" /></a>
quantifies the association between the survival risk and methylation loci M, <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{w}_{1}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{w}_{1}" title="\boldsymbol{w}_{1}" /></a> is a p-vector of fixed effect sizes for clinical covariates.
<p align="center"> 
  <a href="https://www.codecogs.com/eqnedit.php?latex=log(h(t|\boldsymbol{M},\boldsymbol{X})/h_{0}(t))=\sum_{k=1}^{K}M_{k}\gamma&space;_{k}^{TE}&plus;\boldsymbol{Xw}_{1}=\boldsymbol{M\gamma}^{TE}&plus;\boldsymbol{Xw}_{1}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?log(h(t|\boldsymbol{M},\boldsymbol{X})/h_{0}(t))=\sum_{k=1}^{K}M_{k}\gamma&space;_{k}^{TE}&plus;\boldsymbol{Xw}_{1}=\boldsymbol{M\gamma}^{TE}&plus;\boldsymbol{Xw}_{1}" title="log(h(t|\boldsymbol{M},\boldsymbol{X})/h_{0}(t))=\sum_{k=1}^{K}M_{k}\gamma _{k}^{TE}+\boldsymbol{Xw}_{1}=\boldsymbol{M\gamma}^{TE}+\boldsymbol{Xw}_{1}" /></a>
  </p>
  <p align="center">
  <a href="https://www.codecogs.com/eqnedit.php?latex=\gamma&space;_{k}^{TE}\sim&space;N(0,\tau&space;_{1})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\gamma&space;_{k}^{TE}\sim&space;N(0,\tau&space;_{1})" title="\gamma _{k}^{TE}\sim N(0,\tau _{1})" /></a>
    </p>
   τ1 is the variance of <a href="https://www.codecogs.com/eqnedit.php?latex=\gamma&space;_{k}^{TE}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\gamma&space;_{k}^{TE}" title="\gamma&space;_{k}^{TE}" /></a>
   
   
   
   
#### Step 2: Linear mixed-effects model testing for the effect of methylations on gene expression ####
<a href="https://www.codecogs.com/eqnedit.php?latex=\alpha&space;_{k}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha&space;_{k}" title="\alpha&space;_{k}" /></a> quantifies the association between the gene expression G and methylation loci M, <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{w}_{2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{w}_{2}" title="\boldsymbol{w}_{2}" /></a> is a p-vector of fixed effect sizes for clinical covariates.
  <p align="center">
  <a href="https://www.codecogs.com/eqnedit.php?latex=G=\sum_{k=1}^{K}M_{k}\alpha&space;_{k}&plus;\boldsymbol{Xw}_{2}&plus;\varepsilon" target="_blank"><img src="https://latex.codecogs.com/gif.latex?G=\sum_{k=1}^{K}M_{k}\alpha&space;_{k}&plus;\boldsymbol{Xw}_{2}&plus;\varepsilon" title="G=\sum_{k=1}^{K}M_{k}\alpha _{k}+\boldsymbol{Xw}_{2}+\varepsilon" /></a>
</p>
 <p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\alpha&space;_{k}\sim&space;N(0,\tau&space;_{2}),&space;\varepsilon&space;\sim&space;N(0,\sigma&space;_{\varepsilon}^{2})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha&space;_{k}\sim&space;N(0,\tau&space;_{2}),&space;\varepsilon&space;\sim&space;N(0,\sigma&space;_{\varepsilon}^{2})" title="\alpha _{k}\sim N(0,\tau _{2}), \varepsilon \sim N(0,\sigma _{\varepsilon}^{2})" /></a>
  </p>
     τ2 is the variance of <a href="https://www.codecogs.com/eqnedit.php?latex=\alpha&space;_{k}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha&space;_{k}" title="\alpha&space;_{k}" /></a>
  
  
  
  
#### Step 3: Cox linear mixed-effects model testing for the effect of gene expression on the survival outcom ####
<a href="https://www.codecogs.com/eqnedit.php?latex=\gamma&space;_{k}^{DE}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\gamma&space;_{k}^{DE}" title="\gamma&space;_{k}^{DE}" /></a> quantifies the association of methylation loci M on survival risk that is mediated via gene expression G, <a href="https://www.codecogs.com/eqnedit.php?latex=\beta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta" title="\beta" /></a> is the effect eize of gene expression G,<a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{w}_{3}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{w}_{3}" title="\boldsymbol{w}_{3}" /></a> is a p-vector of fixed effect sizes for clinical covariates.
  <p align="center">
  <a href="https://www.codecogs.com/eqnedit.php?latex=\small&space;log(h(t|\boldsymbol{M},G,\boldsymbol{X})/h_{0}(t))=\sum_{k=1}^{K}M_{k}\gamma&space;_{k}^{DE}&plus;G\beta&space;&plus;\boldsymbol{Xw}_{3}=\boldsymbol{M\gamma}^{DE}&plus;G\beta&plus;\boldsymbol{Xw}_{3}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\small&space;log(h(t|\boldsymbol{M},G,\boldsymbol{X})/h_{0}(t))=\sum_{k=1}^{K}M_{k}\gamma&space;_{k}^{DE}&plus;G\beta&space;&plus;\boldsymbol{Xw}_{3}=\boldsymbol{M\gamma}^{DE}&plus;G\beta&plus;\boldsymbol{Xw}_{3}" title="\small log(h(t|\boldsymbol{M},G,\boldsymbol{X})/h_{0}(t))=\sum_{k=1}^{K}M_{k}\gamma&space;_{k}^{DE}+G\beta +\boldsymbol{Xw}_{3}=\boldsymbol{M\gamma}^{DE}+G\beta+\boldsymbol{Xw}_{3}" /></a>
  </p>
   <p align="center">
  <a href="https://www.codecogs.com/eqnedit.php?latex=\gamma&space;_{k}^{DE}\sim&space;N(0,\tau&space;_{3})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\gamma&space;_{k}^{DE}\sim&space;N(0,\tau&space;_{3})" title="\gamma _{k}^{DE}\sim N(0,\tau _{3})" /></a>
    </p>
       τ3 is the variance of <a href="https://www.codecogs.com/eqnedit.php?latex=\gamma&space;_{k}^{DE}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\gamma&space;_{k}^{DE}" title="\gamma&space;_{k}^{DE}" /></a>
      

#### Step 4: Intersection-union survival mixture-adjusted mediation test (IUSMMT) ####
 
IUSMMT verify whether a given gene has mediation effect on the path from methylation CpG sites to the survival risk by testing for:
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=H_{0}:\boldsymbol{\alpha&space;}=\beta&space;=0\Leftrightarrow&space;H_{0}:\tau&space;_{2}=\beta&space;=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?H_{0}:\boldsymbol{\alpha&space;}=\beta&space;=0\Leftrightarrow&space;H_{0}:\tau&space;_{2}=\beta&space;=0" title="H_{0}:\boldsymbol{\alpha }=\beta =0\Leftrightarrow H_{0}:\tau _{2}=\beta =0" /></a>
</p>
This is a joint test including both fixed effect and random effects: the first component of H0 examines the influence of methylation on the gene expression; while the second component examines the impact of gene expression on the survival outcome. Briefly, we derive the test statistic for θ under H0: θ = 0 and τ = 0 as usual, while we derive the score statistic for τ under τ = 0 but without the constraint of θ = 0. By doing this, we ensure that these two statistics are independent. This strategy substantially eases the development of test statistics for the joint test. In conclusion, under this framework two asymptotically independent statistics can be derived: one for the variance component (i.e. τ2) in the generalized linear model and the other for the fixed effect (i.e. <a href="https://www.codecogs.com/eqnedit.php?latex=\beta" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta" title="\beta" /></a>) in the KM Cox model.

# References
Lin X, Cai T, Wu MC, Zhou Q, Liu G, et al. (2011) Kernel machine SNP-set analysis for censored survival outcomes in genome-wide association studies. Genet Epidemiol 35: 620-631.[DOI: 10.1002/gepi.20610](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.20610)

Wu Michael C, Lee S, Cai T, Li Y, Boehnke M, et al. (2011) Rare-Variant Association Testing for Sequencing Data with the Sequence Kernel Association Test. Am J Hum Genet 89: 82-93.[DOI: 10.1016/j.ajhg.2011.05.029.](https://linkinghub.elsevier.com/retrieve/pii/S0002929711002229)

Therneau TM, Grambsch PM, Pankratz VS (2003) Penalized survival models and frailty. Journal of computational and graphical statistics 12: 156-175.[DOI: 10.1198/1061860031365](https://www.tandfonline.com/doi/abs/10.1198/1061860031365)

Langaas M, Lindqvist BH, Ferkingstad E (2005) Estimating the proportion of true null hypotheses, with application to DNA microarray data. J R Stat Soc Ser B 67: 555-572.[DOI: 10.1016/j.jspi.2020.04.011](https://www.sciencedirect.com/science/article/abs/pii/S0378375820300495)

Zhang H, Zheng Y, Zhang Z, Gao T, Joyce B, et al. (2016) Estimating and testing high-dimensional mediation effects in epigenetic studies. Bioinformatics 32: 3150-3154.[DOI: 10.1093/bioinformatics/btw351](https://academic.oup.com/bioinformatics/article/32/20/3150/2196468)


# Cite
Zhonghe Shao<sup>$</sup>, Ting Wang<sup>$</sup>, Meng Zhang, Zhou Jiang, Shuiping Huang<sup>#</sup> and Ping Zeng<sup>#</sup> (2021). IUSMMT: survival mediation analysis of gene expression with multiple DNA methylation exposures and its application to cancers of TCGA.


# Contact
We are very grateful to any questions, comments, or bugs reports; and please contact [Ping Zeng](https://github.com/biostatpzeng) via zpstat@xzhmu.edu.cn.


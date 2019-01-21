
nagelkerke = 
  function(fit, null=NULL, restrictNobs=FALSE)
  {
    TOGGLE =   (class(fit)[1]=="lm"
                | class(fit)[1]=="gls"
                | class(fit)[1]=="lme"
                | class(fit)[1]=="glm"
                | class(fit)[1]=="negbin"
                | class(fit)[1]=="zeroinfl"
                | class(fit)[1]=="clm"
                | class(fit)[1]=="vglm"
                | class(fit)[1]=="betareg"
                | class(fit)[1]=="rq")
    BOGGLE =   (class(fit)[1]=="nls"
                | class(fit)[1]=="lmerMod"
                | class(fit)[1]=="glmerMod"
                | class(fit)[1]=="merModLmerTest"
                | class(fit)[1]=="clmm")
    SMOGGLE =   (class(fit)[1]=="lmerMod"
                 | class(fit)[1]=="glmerMod"
                 | class(fit)[1]=="merModLmerTest"
                 | class(fit)[1]=="vglm")
    ZOGGLE  = (class(fit)[1]=="zeroinfl")
    ZOGGLE2 = (class(fit)[1]=="rq")
    NOGGLE = is.null(null)
    ERROR  = "Note: For models fit with REML, these statistics are based on refitting with ML"
    ERROR2 = "None"
    
    if(!restrictNobs & NOGGLE  & TOGGLE){null = update(fit, ~ 1)}
    if(restrictNobs  & NOGGLE  & TOGGLE){null = update(fit, ~ 1, data=fit$model)}
    
    if(restrictNobs  & !NOGGLE){null = update(null, data=fit$model)}
    
    if(NOGGLE & BOGGLE)
    {ERROR = "You need to supply a null model for nls, lmer, glmer, or clmm"}
    if((!TOGGLE) & (!BOGGLE))
    {ERROR = "This function will work with lm, gls, lme, lmer, glmer, glm, negbin, zeroinfl, nls, clm, clmm, and vglm"}
    
    SMOGGLE2 = (class(null)[1]=="lmerMod"
                | class(null)[1]=="glmerMod"
                | class(null)[1]=="merModLmerTest"
                | class(null)[1]=="vglm")   
    
    Y = matrix(rep(NA,2),
               ncol=1)
    colnames(Y) = ""
    rownames(Y) = c("Model:", "Null:")
    
    Z = matrix(rep(NA, 3),
               ncol=1)
    colnames(Z) = c("Pseudo.R.squared")
    rownames(Z) = c("McFadden", "Cox and Snell (ML)", 
                    "Nagelkerke (Cragg and Uhler)")
    
    X = matrix(rep(NA,4),
               ncol=4)
    colnames(X) = c("Df.diff","LogLik.diff","Chisq","p.value")
    rownames(X) = ""
    
    U = matrix(rep(NA,2),
               ncol=1)
    colnames(U) = ""
    rownames(U) = c("Model:", "Null:")
    
    if(TOGGLE | BOGGLE){
      if (!SMOGGLE){Y[1]= toString(fit$call)}
      if (SMOGGLE){Y[1]= toString(fit@call)}
    }
    
    if(TOGGLE | (BOGGLE & !NOGGLE)){
      
      if (!SMOGGLE2){Y[2]= toString(null$call)}
      if (SMOGGLE2){Y[2]= toString(null@call)}
      
      if(!ZOGGLE & !ZOGGLE2){N = nobs(fit)
                             U[1,1]= nobs(fit); U[2,1]= nobs(null)}
      if(!ZOGGLE &  ZOGGLE2){N = length(fit$y)
                             U[1,1]= length(fit$y); U[2,1]= length(null$y)}
      if(ZOGGLE){N = fit$n
                 U[1,1]= fit$n; U[2,1]= null$n}
      
      if(U[1,1] != U[2,1]){
        ERROR2 = "WARNING: Fitted and null models have different numbers of observations"}
      
      m = suppressWarnings(logLik(fit, REML=FALSE))[1]
      n = suppressWarnings(logLik(null, REML=FALSE))[1]
      mf = 1 - m/n
      Z[1,] = signif(mf, digits=6)
      cs = 1 - exp(-2/N * (m - n))
      Z[2,] = signif(cs, digits=6)
      nk = cs/(1 - exp(2/N * n))
      Z[3,] = signif(nk, digits=6)
      
      o = n - m
      dfm = attr(logLik(fit),"df")
      dfn = attr(logLik(null),"df")
      if(class(fit)[1]=="vglm"){dfm=df.residual(fit)}
      if(class(fit)[1]=="vglm"){dfn=df.residual(null)}
      dff = dfn - dfm
      CHI = 2 * (m - n)
      P = pchisq(CHI, abs(dff), lower.tail = FALSE)
      
      X [1,1] = dff
      X [1,2] = signif(o, digits=5)             
      X [1,3] = signif(CHI, digits=5)
      X [1,4] = signif(P, digits=5)     
    }
    
    W=ERROR
    
    WW=ERROR2
    
    V = list(Y, Z, X, U, W, WW) 
    names(V) = c("Models", "Pseudo.R.squared.for.model.vs.null", 
                 "Likelihood.ratio.test", "Number.of.observations",
                 "Messages", "Warnings")
    return(V)            
  }
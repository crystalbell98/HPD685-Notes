# Linear Longitudinal Modeling Notes

## 0. Unified Notation

### 0.1 Basic Indices and Variables

- $i$: individual
- $t,s$: time points
- $y_{it}$: observed outcome for individual $i$ at time point $t$
- $\mathbf y_i = (y_{i1}, \dots, y_{iT_i})'$: full longitudinal observation vector for individual $i$
- $x_{it}$: time score / time metric (e.g., $0,1,2,3$; can also be age, days since treatment, etc.)
- $w_i$: time-invariant covariate (TIC)
- $z_{it}$: time-varying covariate (TVC)
- $\Delta y_{it} = y_{it} - y_{i,t-1}$: change from $t-1$ to $t$

### 0.2 Unified Linear Growth Notation

Throughout this document, "linear growth" consistently refers to the following equation:

$$
y_{it} = \mu_0 + \mu_1 x_{it} + b_{0i} + b_{1i}x_{it} + \varepsilon_{it}
$$

where:

- $\mu_0$: population mean intercept
- $\mu_1$: population mean slope
- $b_{0i}$: random intercept deviation for individual $i$
- $b_{1i}$: random slope deviation for individual $i$
- $\varepsilon_{it}$: residual / measurement error at time point $t$

The random effects vector is denoted:

$$
\mathbf b_i =
\begin{bmatrix}
b_{0i}\\
b_{1i}
\end{bmatrix},
\qquad
\mathbf D = \text{Var}(\mathbf b_i) =
\begin{bmatrix}
\tau_{00} & \tau_{01}\\
\tau_{01} & \tau_{11}
\end{bmatrix}
$$

- $\tau_{00}$: random intercept variance
- $\tau_{11}$: random slope variance
- $\tau_{01}$: covariance between the random intercept and random slope

### 0.3 Correspondence with SEM Notation

In the longitudinal SEM / LGM literature, the covariance matrix of the random growth factors is commonly denoted $\Psi$, and the residual covariance matrix is denoted $\Theta$. Throughout this document, the correspondence is understood as:

- Random effects covariance matrix in MLM/LMM: $\mathbf D$
- Latent growth factor covariance matrix in SEM/LGM: $\Psi$

Under a simple linear growth model, the two refer to the same quantity and differ only in notation:

$$
\mathbf D \equiv \Psi
$$

Similarly, the Level-1 residual structure $\mathbf R_i$ in MLM/LMM corresponds to the indicator residual structure $\Theta_i$ in SEM.

---

## 1. Methodological Rationale for Longitudinal Research: Why Longitudinal Data Cannot Be Analyzed Directly with OLS

### 1.1 Violation of the Independence Assumption

The core prerequisite of classical OLS is that observations are independently and identically distributed (i.i.d.), such that the residual covariance matrix satisfies:

$$
\Omega = \sigma^2 I
$$

Longitudinal data violate this assumption by design: observations from the same individual at different time points are influenced by the same genetic, personality, environmental, and historical factors, and thus exhibit systematic dependence.

### 1.2 Consequences of Applying OLS Naively

- **Positive autocorrelation**: standard errors are systematically underestimated, inflating $t$ / $F$ statistics and increasing Type I error rates
- **Heterogeneity misattributed to noise**: genuine stable differences between individuals are incorrectly absorbed into random error
- **Confounding of within-person change and between-person differences**: it becomes impossible to distinguish "who started higher" from "who grew faster"
- **Poor handling of missing data**: OLS is ill-suited for unbalanced data structures and is prone to sample bias and information loss

### 1.3 What Longitudinal Modeling Actually Addresses

Longitudinal models do not eliminate dependence; they **parameterize** it:
the failure of the independence assumption is transformed into a question of how the covariance matrix is modeled.

---

## 2. Classical Approaches to Two-Wave Change: Change Score vs. Residualized Change

This section integrates content from the original notes on **raw change scores, residualized change scores, Lord's Paradox, and the measurement and causal limitations of change scores**.

### 2.1 Raw Change Score (Gain Score)

The most direct definition of change is:

$$
\Delta y_i = y_{i2} - y_{i1}
$$

If the change score is used as the outcome variable and regressed on a baseline predictor $w_i$, the model is:

$$
\Delta y_i = \beta_0 + \beta_1 w_i + e_i
$$

**The question this answers:**

- Whose outcome changed the most?
- Which individuals experienced the largest **absolute change**?

**Advantages:**

- Intuitive
- Directly reflects individual-level pre–post change in randomized controlled trials (RCTs)

**Limitations:**

- Influenced by baseline level
- When two measurements are highly correlated, the reliability of the change score may actually be lower than that of the raw scores ("reliability paradox")

### 2.2 Residualized Change Score (ANCOVA)

An alternative approach is to predict the posttest from the pretest, and treat the portion "unexplained by the initial state" as change.

The model is written as:

$$
y_{i2} = \beta_0 + \rho y_{i1} + \beta_1 w_i + e_i
$$

If one is interested only in "deviation from a given baseline," one may first run the regression:

$$
\hat y_{i2} = \hat\beta_0 + \hat\rho y_{i1},
\qquad
r_i = y_{i2} - \hat y_{i2}
$$

Here, $r_i$ is the residualized change.

**The question this answers:**

- Among individuals with the same baseline level, whose posttest performance **exceeded the average expectation**?
- Whose change represents a "deviation from expectation" rather than simply "how much they changed in absolute terms"?

### 2.3 The Essential Distinction Between the Two Change Scores

The two approaches can be summarized in a single sentence:

- **change score** = **absolute change**
- **residualized change** = **baseline-adjusted deviation from expected posttest**

That is:

- $\Delta y_i = 0$ means "no absolute change"
- $r_i = 0$ does not mean "no change"; it means "this person's change is exactly equal to the average expected change for individuals with the same baseline"

### 2.4 A Concrete Example: Why Two Students with Equal Raw Gains Have Different Residualized Changes

Suppose the class-level pretest–posttest regression line is:

$$
\hat y_{\text{post}} = 20 + 0.7\, y_{\text{pre}}
$$

The slope $0.7 < 1$ reflects regression to the mean.

- **Student A**: pre = 40, post = 50; predicted = 48; raw change = $+10$; residualized change = $+2$
- **Student B**: pre = 80, post = 90; predicted = 76; raw change = $+10$; residualized change = $+14$
- **Student C**: pre = 70, post = 69
  predicted = 69
  raw change = $-1$
  residualized change = $0$

Interpretation:

- A and B show identical absolute gains of $+10$
- But B improved far more relative to the average expectation for individuals with the same baseline
- C declined by $1$ point in absolute terms, but this is exactly the average change expected for individuals with the same baseline; hence residualized change = 0

### 2.5 Lord's Paradox: Two Methods, Opposite Conclusions from the Same Data

Lord (1967) noted that, applied to the same dataset, the change score approach and ANCOVA can yield diametrically opposite conclusions.

The classic example involves the effect of a university dining hall's meal plan on the weights of male and female students. At both baseline and follow-up, the overall weight distributions of males and females appear identical.

- **Change score conclusion**: average change is zero for both sexes, therefore "no sex difference"
- **ANCOVA conclusion**: after controlling for baseline, males with the same initial weight are expected to have a higher posttest weight than females, therefore "the meal plan has a significant positive effect for males"

**The underlying reason:**

- The change score approach implicitly fixes the coefficient on the pretest at 1
- ANCOVA allows this coefficient $\rho$ to be estimated from the data (typically $< 1$), and thus "adjusts" for initial differences

In non-randomized observational studies, the two approaches do not target the same estimand, so it is not surprising that they conflict.

### 2.6 When to Use Which Approach?

#### Randomized Controlled Trials (RCTs)

The two approaches tend to converge:

- Random assignment makes baseline values statistically equivalent across groups
- Regression to the mean cancels out across groups
- ANCOVA typically achieves greater statistical power

#### Observational Studies (Non-random Group Assignment)

Caution is required:

- If you are interested in **total effects / absolute change** → change score is more appropriate
- If you are interested in **relative change after accounting for starting point** → ANCOVA / residualized change is more appropriate

However, the prerequisite is: you must believe that "controlling for baseline" is consistent with your causal structure; otherwise, regression-to-the-mean (RTM) bias is easily introduced.

### 2.7 Why Change Scores Can Be Misleading: Measurement-Level Reasons

Suppose the two measurements are:

$$
y_{i1} = T_{i1} + e_{i1},
\qquad
y_{i2} = T_{i2} + e_{i2}
$$

Then the change score is:

$$
\Delta y_i = y_{i2} - y_{i1}
= (T_{i2} - T_{i1}) + (e_{i2} - e_{i1})
$$

If the two errors are independent, then:

$$
\text{Var}(e_{i2} - e_{i1}) = \text{Var}(e_{i2}) + \text{Var}(e_{i1})
$$

That is, measurement errors in the difference **accumulate directly**, so change scores are typically less stable than single-occasion measurements.

In addition, change scores exhibit a mechanical negative correlation with the baseline:

$$
\text{Cov}(\Delta y_i, y_{i1})
=
\text{Cov}(y_{i2} - y_{i1}, y_{i1})
=
\text{Cov}(y_{i2}, y_{i1}) - \text{Var}(y_{i1})
$$

As long as the cross-time correlation is less than 1, this covariance tends to be negative:

- Individuals with a higher baseline are more likely to appear to have gained less
- Individuals with a lower baseline are more likely to appear to have gained more

This does not necessarily reflect genuine psychological or behavioral mechanisms.

### 2.8 Why Change Scores Can Be Misleading: Causal Reasons

In non-randomized data, regressing a change score as the outcome on some baseline exposure often does not correspond to a well-defined causal effect.

This is particularly problematic when the baseline outcome simultaneously serves as a:

- confounder, or
- mediator

In such cases, the regression coefficient for the change score may even be in the opposite direction from the true total effect / direct effect.

**Conclusion**: change scores are not inherently unusable, but in non-randomized causal problems they frequently do not correspond to the effect one intends to estimate.

### 2.9 What Questions Are These Methods Best Suited to Answer?

- "How much did it change?" → both change score and ANCOVA are appropriate
- "What is the shape of the change trajectory?" → two-wave methods are insufficient
- "Do different individuals change at different rates?" → MLM / LGM is required
- "How does a prior state drive subsequent change?" → CLPM or LCSM is required

This is also why the distinction between "absolute change vs. conditional change" will be revisited in the subsequent discussions of CLPM, LGM, and LCSM.

---

## 3. A Unified Mathematical Perspective: Explicitly Modeling Dependence

This section integrates content from the original notes on **the matrix representation of LMM, the implied covariance structure in SEM, and the translation dictionary and equivalence between MLM and SEM/LGM**.

### 3.1 Matrix Form of MLM / LMM

For individual $i$, the LMM is written as:

$$
\mathbf y_i = X_i\beta + Z_i\mathbf b_i + \boldsymbol\varepsilon_i
$$

where:

- $X_i\beta$: fixed effects, representing the population mean structure
- $Z_i\mathbf b_i$: random effects, representing individual deviations from the mean trajectory
- $\boldsymbol\varepsilon_i$: residuals

The marginal covariance matrix implied by this model is:

$$
V_i = \text{Cov}(\mathbf y_i) = Z_i\mathbf D Z_i' + R_i
$$

- $\mathbf D$: random effects covariance matrix, corresponding to heterogeneity in intercepts and slopes
- $R_i$: within-individual residual structure, representing measurement error or short-term fluctuation

This is how LMM "addresses dependence" at the matrix level:
rather than requiring dependence to disappear, it explicitly specifies its sources.

### 3.2 The Same Idea at the Scalar Level

For the linear growth model:

$$
y_{it} = \mu_0 + \mu_1x_{it} + b_{0i} + b_{1i}x_{it} + \varepsilon_{it}
$$

the covariance between two time points $t$ and $s$ for the same individual is:

$$
\text{Cov}(y_{it}, y_{is})
=
\tau_{00}
+ x_{it}\tau_{01}
+ x_{is}\tau_{01}
+ x_{it}x_{is}\tau_{11}
+ \text{Cov}(\varepsilon_{it}, \varepsilon_{is})
$$

If Level-1 residuals are assumed to be independent conditional on the random effects, then:

$$
\text{Cov}(\varepsilon_{it}, \varepsilon_{is}) = 0
\qquad (t \neq s)
$$

and the dependence arises entirely from the shared random intercept and random slope.

In the simplest random-intercept model:

$$
\text{Cov}(y_{it}, y_{is}) = \tau_{00}
\qquad (t \neq s)
$$

Under intensive longitudinal measurement, it is often necessary to additionally model residual serial dependence, for example AR(1):

$$
\text{Cov}(\varepsilon_{it}, \varepsilon_{is})
=
\sigma^2 \rho^{|x_{it} - x_{is}|}
$$

This indicates that errors at nearby time points are more similar, with the correlation weakening as the time interval increases.

### 3.3 Matrix Form of Longitudinal SEM / LGM

In the SEM framework, the longitudinal growth model is written as:

$$
\mathbf y_i = \Lambda_i \eta_i + \varepsilon_i,
\qquad
\eta_i = \alpha + \zeta_i
$$

Therefore:

$$
E(\mathbf y_i) = \Lambda_i \alpha
$$

$$
\Sigma_i = \text{Var}(\mathbf y_i) = \Lambda_i \Psi \Lambda_i' + \Theta_i
$$

where:

- $\Lambda_i$: factor loading matrix, defining the temporal shape (linear, quadratic, latent basis, etc.)
- $\Psi$: latent growth factor covariance matrix
- $\Theta_i$: indicator residual covariance matrix

For the $t$-th time point, if the loading row is $\lambda_t = [1, x_t]$, then the covariance between any two time points is:

$$
\text{Cov}(y_{it}, y_{is})
=
\psi_{00}
+ x_t\psi_{01}
+ x_s\psi_{01}
+ x_tx_s\psi_{11}
+ \theta_{ts}
$$

where:

- $\psi_{00}, \psi_{11}, \psi_{01}$: variances and covariance of the latent intercept and latent slope
- $\theta_{ts}$: residual covariance, used to represent residual dependence of the same indicator across time

### 3.4 From Conditional to Marginal Representation: MLM and SEM Are Mathematically the Same

After integrating out the random effects, the marginal distribution of the LMM is:

$$
E(\mathbf y_i) = X_i\beta
$$

$$
\text{Var}(\mathbf y_i) = Z_i\mathbf D Z_i' + R_i
$$

The marginal structure implied by SEM/LGM is:

$$
E(\mathbf y_i) = \Lambda_i\alpha
$$

$$
\text{Var}(\mathbf y_i) = \Lambda_i\Psi\Lambda_i' + \Theta_i
$$

The two forms are isomorphic.

### 3.5 Translation Dictionary

- $Z_i$ (MLM random-effects design matrix) $\leftrightarrow$ $\Lambda_i$ (SEM factor loading matrix)
- $\mathbf b_i$ (random intercept / slope) $\leftrightarrow$ $\zeta_i$ (individual deviation of the latent growth factor)
- $\mathbf D$ $\leftrightarrow$ $\Psi$
- $R_i$ $\leftrightarrow$ $\Theta_i$
- $X_i\beta$ $\leftrightarrow$ $\Lambda_i\alpha$

The most important thing to remember is:

> **The row vector $[1, x_{it}]$ in MLM corresponds to a row of the factor loading matrix $\Lambda_i$ in SEM.**

### 3.6 Concrete Correspondence for a Four-Wave Linear Growth Model

For four equally spaced time points $t = 0,1,2,3$:

$$
Z_i = \Lambda_i =
\begin{bmatrix}
1 & 0\\
1 & 1\\
1 & 2\\
1 & 3
\end{bmatrix}
$$

- **MLM reading**:

$$
y_{it} = \mu_0 + \mu_1 t + b_{0i} + b_{1i}t + \varepsilon_{it}
$$

- **SEM reading**:

$$
\mathbf y_i = \Lambda_i \eta_i + \varepsilon_i,
\qquad
\eta_i =
\begin{bmatrix}
\mu_0 + b_{0i}\\
\mu_1 + b_{1i}
\end{bmatrix}
$$

As long as the time coding and error structure are consistent, these are two representations of the same growth model.

### 3.7 Random Intercept Only: Full Correspondence with the Null Model

When only a random intercept is included:

$$
\mathbf D = [\tau_{00}]
$$

The corresponding SEM is an intercept-only latent factor with all loadings fixed to 1:

$$
\lambda =
\begin{bmatrix}
1\\
1\\
1
\end{bmatrix}
$$

In this case:

$$
\text{Var}(\mathbf y_i)
=
\lambda \tau_{00}\lambda' + \sigma^2 I
=
\tau_{00}J + \sigma^2 I
$$

This is exactly equivalent to the random-intercept MLM / null model.

### 3.8 Level-2 Predictors in Both Frameworks

Let $w_i$ be a time-invariant covariate.

In MLM, it can predict both the intercept and the slope:

$$
y_{it}
=
(\mu_0 + \beta_0 w_i)
+
(\mu_1 + \beta_1 w_i)x_{it}
+
b_{0i}
+
b_{1i}x_{it}
+
\varepsilon_{it}
$$

In SEM, the same effect is written as:

$$
\eta_{0i} = \mu_0 + \beta_0 w_i + b_{0i}
$$

$$
\eta_{1i} = \mu_1 + \beta_1 w_i + b_{1i}
$$

That is:

- **Cross-level interaction in MLM**
- **"Predicting the latent slope factor" in SEM**

are mathematically the same thing.

### 3.9 When Does the Equivalence Break Down?

This equivalence holds under certain prerequisites:

- Linear Gaussian model
- Consistent time coding
- Consistent residual structure specification
- No additional explicit measurement model

Once any of the following conditions arise, the two approaches diverge:

1. **Different residual structures**

   - LGM more readily allows residual variances to differ across time points by default
   - MLM commonly assumes homogeneous Level-1 residuals by default
2. **SEM explicitly incorporates a measurement model**
   If each time point itself is a latent construct (e.g., depression measured by 5 items), SEM can first specify a measurement model and then a growth model; if MLM still uses sum scores as the outcome, the two are no longer isomorphic.

Therefore, a more precise statement is:

> **The equivalence between MLM and LGM holds for observed-variable simple linear growth models. When SEM further explicitly models measurement error and measurement structure, it goes beyond the scope of a simple one-to-one translation.**

---

## 4. Three Major Frameworks: GEE, MLM/LMM, and Longitudinal SEM

### 4.1 GEE: A Marginal, Population-Averaged Perspective

GEE focuses on:

- How the **population-averaged response** changes when a covariate changes by one unit

It approximates within-cluster dependence by specifying a "working correlation matrix" $R(\alpha)$, and uses a **sandwich variance estimator** to correct the standard errors.

**Advantages:**

- Even if the working correlation matrix is misspecified, the parameter estimates remain consistent as long as the mean model is correctly specified
- Particularly well-suited for public health and policy evaluation contexts where the interest lies in population-level average trends

**Limitations:**

- Less effective at characterizing sources of individual variation
- Missing data typically requires MCAR; handling missingness under MAR often requires multiple imputation (MI)

### 4.2 MLM / LMM: A Conditional Individual Trajectory Perspective

MLM / LMM focuses on:

- How a given individual's trajectory changes, conditional on that individual's random effects
- How much individuals differ in their starting points and rates of change

**In nonlinear models, the coefficients from LMM and GEE typically differ:**

For example, in logistic regression, the fixed-effect estimates from LMM are generally larger than the marginal effects from GEE, because:

- LMM more closely approximates a "pure effect after controlling for individual heterogeneity"
- The GEE effect is attenuated by between-individual heterogeneity

### 4.3 Longitudinal SEM: A Latent Process Perspective

The key advantage of longitudinal SEM is that it:

- Can define a latent construct using multiple indicators
- Explicitly separates "true change" from "measurement error" in growth modeling
- Readily accommodates complex pathways involving mediation, moderation, and parallel process growth

This is its principal advantage over LMM alone.

### 4.4 One-Page Comparison

- **Level of interpretation**

  - GEE = population-averaged trend
  - MLM/LMM = individual conditional trajectory
  - Longitudinal SEM = latent process trajectory
- **Data format**

  - GEE / MLM = long format is typically preferred
  - Longitudinal SEM / LGM = wide format is typically preferred
- **Statistical basis**

  - GEE = quasi-likelihood / marginal model
  - MLM = maximum likelihood (ML / REML)
  - Longitudinal SEM = covariance structure analysis
- **Missing data mechanism**

  - GEE = relies more heavily on MCAR (or requires MI)
  - MLM = can be estimated directly under MAR
  - Longitudinal SEM = can use FIML under MAR
- **Measurement error**

  - GEE = not explicitly addressed
  - MLM = treated as part of the residual
  - Longitudinal SEM = can be explicitly modeled to "purify" the latent variable

---

## 5. MLM / LMM: Conditional Model for Individual Trajectories

### 5.1 Basic Framework: Time at Level 1, Individuals at Level 2

The fundamental idea of longitudinal MLM is:

- **Level 1**: within-individual observations varying over time
- **Level 2**: between-individual differences

The most basic linear growth model is:

$$
y_{it} = \mu_0 + \mu_1 x_{it} + b_{0i} + b_{1i}x_{it} + \varepsilon_{it}
$$

This can also be understood as:

- Each person has their own starting point: $\mu_0 + b_{0i}$
- Each person also has their own rate of growth: $\mu_1 + b_{1i}$

### 5.2 Null Model / Random-Effects ANOVA: First Assessing the Strength of Dependence

The simplest unconditional model is:

$$
y_{it} = \mu_0 + b_{0i} + \varepsilon_{it}
$$

where:

$$
b_{0i} \sim N(0, \tau_{00}),
\qquad
\varepsilon_{it} \sim N(0, \sigma^2)
$$

#### Three Core Implications

1. **Expected value**

$$
E(y_{it}) = \mu_0
$$

1. **Total variance**

$$
\text{Var}(y_{it}) = \tau_{00} + \sigma^2
$$

1. **Covariance between two measurements from the same individual**

$$
\text{Cov}(y_{it}, y_{is}) = \tau_{00}
\qquad (t \neq s)
$$

As long as $\tau_{00} > 0$, repeated measurements from the same individual are not independent.

#### Intraclass Correlation Coefficient (ICC)

$$
\text{ICC} = \frac{\tau_{00}}{\tau_{00} + \sigma^2}
$$

The ICC indicates what proportion of total variance is attributable to between-individual differences, and simultaneously quantifies the strength of data dependence.

- Higher ICC means repeated measurements are less independent
- Higher ICC means OLS is less appropriate

#### Explicit Covariance Matrix (for 3 Time Points)

$$
V_i =
\begin{bmatrix}
\tau_{00}+\sigma^2 & \tau_{00} & \tau_{00}\\
\tau_{00} & \tau_{00}+\sigma^2 & \tau_{00}\\
\tau_{00} & \tau_{00} & \tau_{00}+\sigma^2
\end{bmatrix}
$$

- Diagonal entries: total variance
- Off-diagonal entries: covariance between two measurements from the same individual

This is why OLS fails:
OLS implicitly assumes $V_i = \sigma^2 I$, i.e., all off-diagonal entries are zero.

#### Why Is It Also Called "Random-Effects ANOVA"?

In general clustered data, this model can also be written as a "random-effects ANOVA":

- Traditional ANOVA treats group membership as a fixed factor
- Random-effects ANOVA treats group membership as a set of clusters randomly sampled from a population

In the longitudinal setting, the cluster is the "individual," so the same logic applies directly.

#### Limitations of the ICC

The ICC is not reliable in all situations:

- **Binary outcomes**: in models such as logistic regression, there is no simple, unified Level-1 variance term, so the ICC is no longer straightforwardly defined
- **Models with random slopes**: within-cluster correlation varies as a function of the predictor value, and a single ICC is insufficient as a summary
- **Very small clusters**: population mean estimates are unstable, and ICC estimates tend to be biased

### 5.3 Random Intercept and Random Slope Models

#### Random Intercept Model

If only starting points are allowed to vary while the rate of change is fixed:

$$
y_{it} = \mu_0 + \mu_1 x_{it} + b_{0i} + \varepsilon_{it}
$$

This indicates that each person has a different baseline, but the time effect is fixed.

#### Random Slope Model

If individuals differ not only in their starting points but also in their rates of change:

$$
y_{it} = \mu_0 + \mu_1 x_{it} + b_{0i} + b_{1i}x_{it} + \varepsilon_{it}
$$

In this case:

- $b_{0i}$: accounts for "who starts higher"
- $b_{1i}$: accounts for "who changes faster"

### 5.4 How Are the Parameters Interpreted?

- **$\mu_0$**: population mean level at time zero
- **$\mu_1$**: population mean rate of change
- **$\tau_{00}$**: magnitude of between-individual differences in starting point
- **$\tau_{11}$**: magnitude of between-individual differences in rate of growth
- **$\tau_{01}$**: whether starting point and rate of growth are correlated

The interpretation of $\tau_{01}$ is often of particular research interest:

- **Positive correlation**: individuals who start higher subsequently grow faster (Matthew effect)
- **Negative correlation**: individuals who start higher subsequently grow more slowly (ceiling / plateau effect)

### 5.5 Four Key Quantities to Examine When Reading MLM Results

- **Fixed time effect**: average growth / decline trend
- **Random intercept variance**: between-individual differences in starting point
- **Random slope variance**: between-individual differences in rate of change
- **Intercept–slope covariance**: whether individuals who start higher subsequently accelerate or decelerate

### 5.6 R / lmer Code Example (GPA Setting)

```r
# Random intercept model: each student has a different starting point, slope is fixed
m1 <- lmer(gpa ~ time + (1 | student), data = longdat)

# Random intercept + random slope + sex interaction: starting point and growth rate both vary by individual
m2 <- lmer(gpa ~ time * sex + (time | student), data = longdat)
```

Interpretation:

- `(1 | student)`: allows only the intercept to vary across individuals
- `(time | student)`: both the intercept and the time slope vary across individuals

In LGM terms, "having a variable predict the intercept and slope" corresponds in MLM to "main effect + time interaction + random intercept / random slope."

### 5.7 Applied Example: Adolescent Alcohol Use Trajectories

Alcohol use frequency was measured at ages 14, 15, and 16. After fitting a random slope model:

- **Fixed effects**

  - Intercept = 0.651: mean usage at age 14 is significantly different from zero
  - Time slope = 0.271: average increase of 0.271 units per year
- **Random effects**

  - Intercept variance = 0.624: significant between-individual differences at baseline
  - Slope variance = 0.151: individual differences in rate of growth are present
  - Intercept–slope covariance = −0.07: higher initial use at age 14 is associated with slower subsequent growth
- **Level-2 predictors**

  - COA (children of alcoholics) significantly predicts higher initial use ($\beta = 0.743$)
  - But does not significantly predict the rate of growth

### 5.8 Within-Person / Between-Person Decomposition: Why Separation Is Necessary

Understanding the WP/BP decomposition is central to MLM and longitudinal analysis. Without this decomposition, the estimated "total effect" conflates two fundamentally different processes, potentially resulting in an **Ecological Fallacy**—incorrectly inferring individual-level mechanisms from group-level associations.

- **Between-Person (BP) differences**: stable differences between individuals (trait-level variation)
- **Within-Person (WP) differences**: fluctuations within the same individual across time (state-level change)

For a time-varying covariate $z_{it}$, its coefficient typically conflates:

- **Within-person effect**: how does $y$ change when a person's $z$ is higher than their own usual level?
- **Between-person effect**: do individuals with a consistently higher average $z$ also tend to have higher $y$ on average?

These two effects may operate in the same direction or in opposite directions. A classic example is stress and negative affect:

- **Between-person**: individuals under chronically high stress have higher average negative affect
- **Within-person**: on days when a person's stress exceeds their own mean, is their negative affect also elevated?

#### Person-Mean Centering (Within-Person Centering)

The time-varying covariate $z_{it}$ is decomposed into two **orthogonal** components:

$$
z_{it} = \underbrace{\bar z_i}_{\text{BP}} + \underbrace{(z_{it} - \bar z_i)}_{\text{WP}}
$$

- $\bar z_i$: individual mean, entered at Level 2 to capture **Between-Person** variation (the "who" effect)
- $z_{it} - \bar z_i$: person-mean-centered value, entered at Level 1 to capture **Within-Person** variation (the "when" effect)

Both components are entered simultaneously into MLM to estimate the WP and BP effects separately:

$$
y_{ij} = \gamma_{00} + \underbrace{\gamma_{10}(z_{ij} - \bar z_j)}_{\text{WP effect}} + \underbrace{\gamma_{01}\bar z_j}_{\text{BP effect}} + u_{0j} + r_{ij}
$$

- **$\gamma_{10}$ (Within-Person Effect)**: how $y$ changes when an individual's current $z$ is above their own mean
- **$\gamma_{01}$ (Between-Person Effect)**: whether individuals with a higher average $z$ also have a higher average $y$
- **$u_{0j}$**: random intercept (between-individual baseline differences); **$r_{ij}$**: residual

**Benefits:**

- The Level-1 $\gamma_{10}$ is a pure within-person dynamic effect, uncontaminated by stable individual differences
- If $\gamma_{10} \neq \gamma_{01}$, the WP and BP processes differ, and failing to decompose them leads to severely biased parameter estimates

#### Specific Application Examples

##### Example 1: Exercise and Heart Attack (Simpson's Paradox)

- **Within-Person**: at the individual level, the instantaneous risk of cardiac events is elevated during vigorous exercise
- **Between-Person**: at the population level, individuals who exercise regularly have a substantially lower baseline cardiac risk than sedentary individuals
- **Conclusion**: without decomposition, one might incorrectly infer that "exercise is unrelated to heart disease" or even "exercise is harmful"

##### Example 2: Self-Esteem and Self-Enhancement

A study of 60 students on 14 personality traits, with "trait importance" person-mean-centered:

- **Within-Person effect**: $\gamma_{10} = 0.37$, meaning that within the same student, traits perceived as more important are associated with greater self-enhancement on that trait
- **Substantive implication**: demonstrates that self-enhancement is a context-sensitive psychological mechanism rather than a fixed trait

##### Example 3: GDP and National Well-Being (Easterlin Paradox)

- **Between-Country (BP)**: wealthier countries typically report higher average well-being than poorer countries
- **Within-Country (WP)**: as a country's GDP grows from year to year, national well-being does not necessarily increase correspondingly (the WP effect may be near zero)
- **Substantive implication**: this is the "Easterlin Paradox"—simply pursuing economic growth does not necessarily improve national well-being

> **In summary**: Between-Person addresses "**who**" scores higher; Within-Person addresses "**when**" one scores higher. In longitudinal analysis, person-mean centering is consistently recommended to obtain clean WP estimates while simultaneously including individual means to capture BP estimates.

---

## 6. Longitudinal SEM / LGM: Latent Processes, Time Coding, and Measurement Models

### 6.1 Basic Structure of LGM: Growth Parameters as Latent Variables

In LGM, repeated measurements are typically arranged in wide-format data, with each time point as a manifest variable. Latent growth factors are identified by fixing factor loadings according to the temporal design.

- **Latent intercept factor**: all time-point loadings fixed to 1
- **Latent slope factor**: loadings fixed according to the time scores, e.g., $0,1,2,3$

The factor loading matrix for a four-wave linear growth model is:

$$
\Lambda =
\begin{bmatrix}
1 & 0\\
1 & 1\\
1 & 2\\
1 & 3
\end{bmatrix}
$$

The corresponding latent growth factors can be written as:

$$
\eta_i =
\begin{bmatrix}
\mu_0 + b_{0i}\\
\mu_1 + b_{1i}
\end{bmatrix}
$$

Therefore:

$$
\mathbf y_i = \Lambda \eta_i + \varepsilon_i
$$

The core parameters of interest in LGM correspond almost one-to-one with those in MLM:

- Intercept mean
- Slope mean
- Intercept variance
- Slope variance
- Intercept–slope covariance
- Residual variance at each time point

#### lavaan Code Example

```r
model <- '
  i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4
  s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4
'
fit <- growth(model, data = widedat)
```

### 6.2 LGM Simultaneously Estimates Mean Structure and Covariance Structure

LGM does not merely estimate "whether the average is changing"; it simultaneously estimates:

- **Mean structure**: population-level average growth trajectory
- **Covariance structure**: distribution of individual differences

This is also why LGM is typically not a saturated model:
you impose a functional form (e.g., linear growth), and the model then asks: can this form account for the data?

If fit indices (e.g., CFI, TLI, RMSEA, SRMR) are poor, the imposed growth form is too simple to explain the true pattern of mean shifts and covariance structure.

### 6.3 Time Coding, Centering, and the Meaning of the Intercept

Time coding is not a technical detail; it determines "what moment the intercept represents."

#### In LGM

- If the slope loading at the first time point is set to 0 → the intercept represents the **baseline level**
- If the loading at some intermediate time point is set to 0 → the intercept represents the **mid-development level**

#### In MLM

The same logic applies:

- wherever the zero point of $x_{it}$ is placed
- that is what $\mu_0$ represents

#### Common Time Metrics

- **Wave number**: simple, but the intercept typically represents only "the first wave"
- **Age**: more developmentally meaningful, but different individuals may have non-overlapping age ranges
- **Days since treatment / time since event**: appropriate for clinical intervention studies
- **Years to death**: appropriate for studies of terminal illness

#### What Changes When the Zero Point Is Shifted?

- The direction and unit of the slope are not affected
- But the following change:
  - Intercept mean
  - Intercept variance
  - Intercept–slope covariance

### 6.4 Latent Basis Model (LBM)

To allow for growth that is not strictly linear, one may fix:

- The loading at the first time point = 0
- The loading at the last time point = 1
- The loadings at intermediate time points to be freely estimated from the data

The freely estimated intermediate loadings then represent:

- The proportion of the total change that has been completed by that time point

This approach is well suited for describing developmental trajectories that are "fast early and slow late" or "slow early and fast late."

### 6.5 Longitudinal Measurement Invariance: Verify the Scale Has Not Shifted Before Fitting a Growth Model

If the construct of interest is latent (e.g., depression, motivation, self-esteem), measurement invariance across time points must be established before comparing change over time.

#### Metric Invariance (Weak Invariance)

- Factor loadings are equal across time points
- Ensures that the unit of the latent variable is consistent

#### Scalar Invariance (Strong Invariance)

- Factor loadings and indicator intercepts are equal across time points
- A necessary prerequisite for comparing latent variable means over time

If scalar invariance does not hold, observed "growth" may reflect only:

- Changes in the items themselves
- Indicator drift
- Shifts in how respondents interpret the scale

rather than genuine change in the latent construct.

## 9. Modeling Decisions: Which Type of Longitudinal Model Should Be Used?

### When MLM Is More Appropriate

- Many time points (e.g., $\ge 10$ occasions)
- Irregular time intervals
- Severely unbalanced data
- Need for multiple levels of nesting (students within classrooms, classrooms within schools)
- Smaller samples (e.g., $N < 100$), where REML tends to be more stable

#### When LGM / Longitudinal SEM Is More Appropriate

- Fewer time points (e.g., 3–8 occasions) with regular intervals
- The construct of interest is latent (depression, motivation, ability, etc.)
- Multiple indicators per time point
- Need to explicitly correct for measurement error
- Need to test longitudinal measurement invariance
- Need to evaluate overall model fit (CFI, RMSEA, etc.)
- Uncertainty about the time scores, with a desire to let the data estimate them (e.g., latent basis model)

## 10. One-Page Summary: The Core Thread of This Document

The entire set of longitudinal modeling notes can be compressed into the following thread:

1. **OLS is inappropriate** because longitudinal data violate the independence assumption
2. **Change can first be described with two-wave methods**, but change scores and residualized change scores do not answer the same question
3. **The core of modern longitudinal modeling** is not to avoid dependence but to write it into the covariance structure
4. **GEE / MLM / Longitudinal SEM** correspond to three perspectives: population-averaged, individual trajectory, and latent process
5. **MLM and LGM** are often mathematically equivalent under simple linear growth, but SEM can further explicitly address measurement error
6. **CLPM / RI-CLPM / LCSM** further address "who drives whom" and "how a prior state determines the next change"
7. The true criterion for model selection is always: **research question + data structure + estimation target**

---

## Appendix: Unified Symbol Reference Table

- $y_{it}$: observed value for individual $i$ at time $t$
- $x_{it}$: time score / time metric
- $w_i$: time-invariant covariate (TIC)
- $z_{it}$: time-varying covariate (TVC)
- $\Delta y_{it}$: change score
- $\mu_0$: population mean intercept
- $\mu_1$: population mean slope
- $b_{0i}$: random intercept deviation
- $b_{1i}$: random slope deviation
- $\tau_{00}$: random intercept variance
- $\tau_{11}$: random slope variance
- $\tau_{01}$: random intercept–random slope covariance
- $\varepsilon_{it}$: residual / measurement error
- $\mathbf D$: random effects covariance matrix (denoted $\Psi$ in the SEM literature)
- $R_i$: residual covariance matrix (denoted $\Theta_i$ in the SEM literature)
- $\Lambda$: factor loading matrix / time score matrix
- $\eta_i$: latent growth factor vector
- $\Psi$: latent growth factor covariance matrix
- $\Theta_i$: indicator residual covariance matrix

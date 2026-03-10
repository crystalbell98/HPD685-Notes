# Linear Longitudinal Modeling Notes

## Unified Notation

### Basic Indices and Variables

- $i$: individual
- $t,s$: time points
- $y_{it}$: observed outcome for individual $i$ at time point $t$
- $\mathbf y_i = (y_{i1}, \dots, y_{iT_i})'$: full longitudinal observation vector for individual $i$
- $x_{it}$: time score / time metric (e.g., $0,1,2,3$; can also be age, days since treatment, etc.)
- $\mathbf D$: random effects covariance matrix (= $\Psi$ ,latent growth factor covariance matrix in the SEM literature)
- $R_i$: residual covariance matrix (= $\Theta_i$ indicator residual covariance matrix in the SEM literature)
- $\Lambda$: factor loading matrix / time score matrix
- $\eta_i$: latent growth factor vector

### Unified Linear Growth Notation

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

### Level-2 Predictors 

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


In the SEM framework, the longitudinal growth model is written as:

$$
\mathbf y_i = \Lambda_i \eta_i + \varepsilon_i,
\qquad
\eta_i = \alpha + \zeta_i
$$


the same level-2 effect is written as:

$$
\eta_{0i} = \mu_0 + \beta_0 w_i + b_{0i}
$$

$$
\eta_{1i} = \mu_1 + \beta_1 w_i + b_{1i}
$$


## Methodological Rationale for Longitudinal Research: Why Longitudinal Data Cannot Be Analyzed Directly with OLS

### Violation of the Independence Assumption

The core prerequisite of classical OLS is that observations are independently and identically distributed (i.i.d.), such that the residual covariance matrix satisfies:

$$
\Omega = \sigma^2 I
$$

Longitudinal data violate this assumption by design: observations from the same individual at different time points are influenced by the same genetic, personality, environmental, and historical factors, and thus exhibit systematic dependence.

### What Longitudinal Modeling Actually Addresses

Longitudinal models do not eliminate dependence; they **parameterize** it:
the failure of the independence assumption is transformed into a question of how the covariance matrix is modeled.

---

## A Unified Mathematical Perspective: Explicitly Modeling Dependence

This section integrates content from the original notes on **the matrix representation of LMM, the implied covariance structure in SEM, and the translation dictionary and equivalence between MLM and SEM/LGM**.

### Matrix Form of MLM / LMM

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

### The Same Idea at the Scalar Level

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

### Matrix Form of Longitudinal SEM / LGM

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

### MLM and SEM Are Mathematically the Same

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

### Translation Dictionary

- $Z_i$ (MLM random-effects design matrix) $\leftrightarrow$ $\Lambda_i$ (SEM factor loading matrix)
- $\mathbf b_i$ (random intercept / slope) $\leftrightarrow$ $\zeta_i$ (individual deviation of the latent growth factor)
- $\mathbf D$ $\leftrightarrow$ $\Psi$
- $R_i$ $\leftrightarrow$ $\Theta_i$
- $X_i\beta$ $\leftrightarrow$ $\Lambda_i\alpha$

The most important thing to remember is:

> **The row vector $[1, x_{it}]$ in MLM corresponds to a row of the factor loading matrix $\Lambda_i$ in SEM.**

### Concrete Correspondence for a Four-Wave Linear Growth Model

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

### Random Intercept Only: Full Correspondence with the Null Model

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


### When Does the Equivalence Break Down?

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

## MLM / LMM: Conditional Model for Individual Trajectories

### Basic Framework: Time at Level 1, Individuals at Level 2

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

### Null Model / Random-Effects ANOVA: First Assessing the Strength of Dependence

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

### Random Intercept and Random Slope Models

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

#### How Are the Parameters Interpreted?

- **$\mu_0$**: population mean level at time zero
- **$\mu_1$**: population mean rate of change
- **$\tau_{00}$**: magnitude of between-individual differences in starting point
- **$\tau_{11}$**: magnitude of between-individual differences in rate of growth
- **$\tau_{01}$**: whether starting point and rate of growth are correlated

The interpretation of $\tau_{01}$ is often of particular research interest:

- **Positive correlation**: individuals who start higher subsequently grow faster (Matthew effect)
- **Negative correlation**: individuals who start higher subsequently grow more slowly (ceiling / plateau effect)

---

## One-Page Summary: The Core Thread of This Document

The entire set of longitudinal modeling notes can be compressed into the following thread:

1. **OLS is inappropriate** because longitudinal data violate the independence assumption
2. **The core of modern longitudinal modeling** is not to avoid dependence but to write it into the covariance structure
3. **GEE / MLM / Longitudinal SEM** correspond to three perspectives: population-averaged, individual trajectory, and latent process
4. **MLM and LGM** are often mathematically equivalent under simple linear growth, but SEM can further explicitly address measurement error

---

## Design Matrix

### Core Concept

Multiple regression is expressed compactly in matrix form as: $Y' = X\hat{\beta}$

- **$Y'$**: a column vector of predicted values
- **$X$ (Design Matrix)**: a matrix containing a column of ones (for the intercept) and columns of observed values for all predictors
- **$\hat{\beta}$**: a vector of estimated regression coefficients

**Expanded form with a single predictor (5 observations):**

$$
Y'_i = \hat{\beta}_0 + \hat{\beta}_1 X_{1i} \;\equiv\; \hat{\beta}_0(1) + \hat{\beta}_1 X_{1i}
$$

$$
\begin{bmatrix} Y'_1 \\ Y'_2 \\ Y'_3 \\ Y'_4 \\ Y'_5 \end{bmatrix} = \begin{bmatrix} 1 & X_{11} \\ 1 & X_{12} \\ 1 & X_{13} \\ 1 & X_{14} \\ 1 & X_{15} \end{bmatrix} \begin{bmatrix} \hat{\beta}_0 \\ \hat{\beta}_1 \end{bmatrix} = \begin{bmatrix} \hat{\beta}_0 + \hat{\beta}_1 X_{11} \\ \hat{\beta}_0 + \hat{\beta}_1 X_{12} \\ \hat{\beta}_0 + \hat{\beta}_1 X_{13} \\ \hat{\beta}_0 + \hat{\beta}_1 X_{14} \\ \hat{\beta}_0 + \hat{\beta}_1 X_{15} \end{bmatrix}
$$

**Expanded form with two predictors:**

$$
Y'_i = \hat{\beta}_0 + \hat{\beta}_1 X_{1i} + \hat{\beta}_2 X_{2i}
$$

$$
\begin{bmatrix} Y'_1 \\ Y'_2 \\ Y'_3 \\ Y'_4 \\ Y'_5 \end{bmatrix} = \begin{bmatrix} 1 & X_{11} & X_{21} \\ 1 & X_{12} & X_{22} \\ 1 & X_{13} & X_{23} \\ 1 & X_{14} & X_{24} \\ 1 & X_{15} & X_{25} \end{bmatrix} \begin{bmatrix} \hat{\beta}_0 \\ \hat{\beta}_1 \\ \hat{\beta}_2 \end{bmatrix} = \begin{bmatrix} \hat{\beta}_0 + \hat{\beta}_1 X_{11} + \hat{\beta}_2 X_{21} \\ \hat{\beta}_0 + \hat{\beta}_1 X_{12} + \hat{\beta}_2 X_{22} \\ \hat{\beta}_0 + \hat{\beta}_1 X_{13} + \hat{\beta}_2 X_{23} \\ \hat{\beta}_0 + \hat{\beta}_1 X_{14} + \hat{\beta}_2 X_{24} \\ \hat{\beta}_0 + \hat{\beta}_1 X_{15} + \hat{\beta}_2 X_{25} \end{bmatrix}
$$

### OLS Matrix Solution

The coefficient vector is obtained in one step via matrix algebra: construct the design matrix $X$, compute its transpose $X'$, form the product $X'X$, invert it, and multiply by $X'Y$:

$$
\hat{\beta} = (X'X)^{-1}X'Y
$$

> As long as $X'X$ is invertible (no perfect multicollinearity), this formula yields the exact solution directly, with no iteration.

### Covariance Matrix of the Coefficient Estimates

The covariance matrix of the vector $\hat{\beta}$:

$$
\text{cov}(\hat{\beta}) = \sigma^2 (X'X)^{-1}
$$

where the residual variance $\sigma^2$ is estimated from the residual vector $\epsilon$:

$$
\sigma^2 = \frac{\epsilon'\epsilon}{n - s}
$$

- **$\epsilon$**: the residual vector ($= Y - X\hat{\beta}$)
- **$s$**: the number of regression coefficients estimated in the model (including the intercept)
- **$n - s$**: the residual degrees of freedom

The diagonal elements of the covariance matrix are the variances of the individual coefficient estimates; their square roots are the **standard errors (SE)** of the respective coefficients. This matrix can be extracted in software via the `vcov()` function to verify the SE values reported in the output.

### Matrix Representation of $\text{Var}(\mathbf{y})$ in OLS

The OLS data-generating model is $\mathbf{y} = X\beta + \boldsymbol{\epsilon}$, where $X\beta$ is fixed (non-random), so:

$$
\text{Var}(\mathbf{y}) = \text{Var}(\boldsymbol{\epsilon}) = \sigma^2 \mathbf{I}_n
$$

Expanded as an $n \times n$ matrix:

$$
\text{Var}(\mathbf{y}) = \begin{bmatrix} \sigma^2 & 0 & \cdots & 0 \\ 0 & \sigma^2 & \cdots & 0 \\ \vdots & & \ddots & \vdots \\ 0 & 0 & \cdots & \sigma^2 \end{bmatrix}
$$

**The diagonal elements are the variance of each observation (all equal); the off-diagonal elements are all zero, indicating mutual independence across observations.**

This is precisely the two core assumptions of OLS expressed in the language of covariance matrices: **Homoscedasticity + Independence**.

#### Correspondence to MLM / LSEM

- **OLS**: $\text{Var}(\mathbf{y}) = \sigma^2 \mathbf{I}_n$ — no random effects; residuals are independent and identically distributed
- **MLM**: $V_i = Z_i \mathbf{D} Z_i' + R_i$ — random effects plus within-individual residual structure
- **LSEM**: $\Sigma(\theta) = \Lambda\Psi\Lambda' + \Theta$ — latent factors plus indicator residuals

OLS is the **most constrained special case** of the three: setting $Z_i = \mathbf{0}$ (no random effects) and $R_i = \sigma^2 \mathbf{I}$ (homogeneous independent residuals) reduces the MLM formula immediately to the OLS form $\sigma^2\mathbf{I}$; likewise, setting $\Lambda = \mathbf{0}$ and $\Theta = \sigma^2\mathbf{I}$ reduces the LSEM to OLS.

#### Extensions When Assumptions Are Violated

- **Heteroscedasticity**: off-diagonal elements remain zero (independence holds), but diagonal elements differ:

$$
\text{Var}(\mathbf{y}) = \begin{bmatrix} \sigma_1^2 & & \\ & \ddots & \\ & & \sigma_n^2 \end{bmatrix} = \Omega
$$

- **Autocorrelation (e.g., AR(1))**: diagonal elements are equal, but off-diagonal elements are non-zero:

$$
\text{Var}(\mathbf{y}) = \sigma^2 \begin{bmatrix} 1 & \rho & \rho^2 & \cdots \\ \rho & 1 & \rho & \cdots \\ \rho^2 & \rho & 1 & \cdots \\ \vdots & & & \ddots \end{bmatrix}
$$

- **General case (GLS)**: $\text{Var}(\mathbf{y}) = \sigma^2\Omega$; the OLS formula becomes the Generalized Least Squares (GLS) estimator:

$$
\hat{\beta}_{\text{GLS}} = (X'\Omega^{-1}X)^{-1}X'\Omega^{-1}Y
$$

When $\Omega = \mathbf{I}$, this reduces back to the OLS estimator $\hat{\beta} = (X'X)^{-1}X'Y$.

---

## Within-Person / Between-Person Decomposition: Why Separation Is Necessary

Understanding the WP/BP decomposition is central to MLM and longitudinal analysis. Without this decomposition, the estimated "total effect" conflates two fundamentally different processes, potentially resulting in an **Ecological Fallacy**—incorrectly inferring individual-level mechanisms from group-level associations.

- **Between-Person (BP) differences**: stable differences between individuals (trait-level variation)
- **Within-Person (WP) differences**: fluctuations within the same individual across time (state-level change)

For a time-varying covariate $z_{it}$, its coefficient typically conflates:

- **Within-person effect**: how does $y$ change when a person's $z$ is higher than their own usual level?
- **Between-person effect**: do individuals with a consistently higher average $z$ also tend to have higher $y$ on average?

These two effects may operate in the same direction or in opposite directions. A classic example is stress and negative affect:

- **Between-person**: individuals under chronically high stress have higher average negative affect
- **Within-person**: on days when a person's stress exceeds their own mean, is their negative affect also elevated?

### Person-Mean Centering (Within-Person Centering)

The time-varying covariate $z_{it}$ is decomposed into two **orthogonal** components:

$$
z_{it} = \underbrace{\bar z_i}_{\text{BP}} + \underbrace{(z_{it} - \bar z_i)}_{\text{WP}}
$$

- $\bar z_i$: individual mean, entered at Level 2 to capture **Between-Person** variation (the "who" effect)
- $z_{it} - \bar z_i$: person-mean-centered value, entered at Level 1 to capture **Within-Person** variation (the "when" effect)

Both components are entered simultaneously into MLM to estimate the WP and BP effects separately:

$$
y_{it} = \mu_0 + \underbrace{\beta_{\text{WP}}(z_{it} - \bar z_i)}_{\text{WP effect}} + \underbrace{\beta_{\text{BP}}\bar z_i}_{\text{BP effect}} + b_{0i} + \varepsilon_{it}
$$

- **$\beta_{\text{WP}}$ (Within-Person Effect)**: how $y$ changes when an individual's current $z$ is above their own mean
- **$\beta_{\text{BP}}$ (Between-Person Effect)**: whether individuals with a higher average $z$ also have a higher average $y$
- **$b_{0i}$**: random intercept (between-individual baseline differences); **$\varepsilon_{it}$**: residual

**Benefits:**

- The Level-1 $\beta_{\text{WP}}$ is a pure within-person dynamic effect, uncontaminated by stable individual differences
- If $\beta_{\text{WP}} \neq \beta_{\text{BP}}$, the WP and BP processes differ, and failing to decompose them leads to severely biased parameter estimates

### Specific Application Examples

#### Example 1: Exercise and Heart Attack (Simpson's Paradox)

- **Within-Person**: at the individual level, the instantaneous risk of cardiac events is elevated during vigorous exercise
- **Between-Person**: at the population level, individuals who exercise regularly have a substantially lower baseline cardiac risk than sedentary individuals
- **Conclusion**: without decomposition, one might incorrectly infer that "exercise is unrelated to heart disease" or even "exercise is harmful"

#### Example 2: Self-Esteem and Self-Enhancement

A study of 60 students on 14 personality traits, with "trait importance" person-mean-centered:

- **Within-Person effect**: $\beta_{\text{WP}} = 0.37$, meaning that within the same student, traits perceived as more important are associated with greater self-enhancement on that trait
- **Substantive implication**: demonstrates that self-enhancement is a context-sensitive psychological mechanism rather than a fixed trait

#### Example 3: GDP and National Well-Being (Easterlin Paradox)

- **Between-Country (BP)**: wealthier countries typically report higher average well-being than poorer countries
- **Within-Country (WP)**: as a country's GDP grows from year to year, national well-being does not necessarily increase correspondingly (the WP effect may be near zero)
- **Substantive implication**: this is the "Easterlin Paradox"—simply pursuing economic growth does not necessarily improve national well-being

> **In summary**: Between-Person addresses "**who**" scores higher; Within-Person addresses "**when**" one scores higher. In longitudinal analysis, person-mean centering is consistently recommended to obtain clean WP estimates while simultaneously including individual means to capture BP estimates.

---

## Classical Approaches to Two-Wave Change: Change Score vs. Residualized Change

This section integrates content from the original notes on **raw change scores, residualized change scores, Lord's Paradox, and the measurement and causal limitations of change scores**.

### Raw Change Score (Gain Score)

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

### Residualized Change Score (ANCOVA)

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

### The Essential Distinction Between the Two Change Scores

The two approaches can be summarized in a single sentence:

- **change score** = **absolute change**
- **residualized change** = **baseline-adjusted deviation from expected posttest**

That is:

- $\Delta y_i = 0$ means "no absolute change"
- $r_i = 0$ does not mean "no change"; it means "this person's change is exactly equal to the average expected change for individuals with the same baseline"

### A Concrete Example: Why Two Students with Equal Raw Gains Have Different Residualized Changes

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

### Lord's Paradox: Two Methods, Opposite Conclusions from the Same Data

Lord (1967) noted that, applied to the same dataset, the change score approach and ANCOVA can yield diametrically opposite conclusions.

The classic example involves the effect of a university dining hall's meal plan on the weights of male and female students. At both baseline and follow-up, the overall weight distributions of males and females appear identical.

- **Change score conclusion**: average change is zero for both sexes, therefore "no sex difference"
- **ANCOVA conclusion**: after controlling for baseline, males with the same initial weight are expected to have a higher posttest weight than females, therefore "the meal plan has a significant positive effect for males"

**The underlying reason:**

- The change score approach implicitly fixes the coefficient on the pretest at 1
- ANCOVA allows this coefficient $\rho$ to be estimated from the data (typically $< 1$), and thus "adjusts" for initial differences

In non-randomized observational studies, the two approaches do not target the same estimand, so it is not surprising that they conflict.

### When to Use Which Approach?

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

### Why Change Scores Can Be Misleading: Measurement-Level Reasons

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

assume that the fluctuations (variances) of the two measurements are the same, i.e., Var($y_{i1}$) = Var($y_{i2}$) = σ². Then:

- Correlation coefficient formula: $ρ = Cov(y_{i1}, y_{i2}) / (σ · σ)$
- Substitute into the formula: $Cov(Δy_i, y_{i1}) = ρ σ² - σ² = σ² (ρ - 1)$
- As long as ρ < 1 (as long as the *two measurements are not completely deterministic* ), then (ρ - 1) is necessarily negative.
- This means that Change Score (change value) and Baseline (initial value) inherently have a negative correlation.
- Individuals with a higher baseline are more likely to appear to have gained less
- Individuals with a lower baseline are more likely to appear to have gained more

This does not necessarily reflect genuine psychological or behavioral mechanisms.

### Why Change Scores Can Be Misleading: Causal Reasons

In non-randomized data, regressing a change score as the outcome on some baseline exposure often does not correspond to a well-defined causal effect.

This is particularly problematic when the baseline outcome simultaneously serves as a:

- confounder, or
- mediator

In such cases, the regression coefficient for the change score may even be in the opposite direction from the true total effect / direct effect.

**Conclusion**: change scores are not inherently unusable, but in non-randomized causal problems they frequently do not correspond to the effect one intends to estimate.

### What Questions Are These Methods Best Suited to Answer?

- "How much did it change?" → both change score and ANCOVA are appropriate
- "What is the shape of the change trajectory?" → two-wave methods are insufficient
- "Do different individuals change at different rates?" → MLM / LGM is required
- "How does a prior state drive subsequent change?" → LGM is required

This is also why the distinction between "absolute change vs. conditional change" will be revisited in the subsequent discussions of CLPM, LGM, and LCSM.

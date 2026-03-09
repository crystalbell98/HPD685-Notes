# Statistical Modeling Notes

---

## 1. Ordinary Least Squares Regression (OLS Regression)

### 1.1 Principles

OLS regression examines the linear relationship between a continuous dependent variable $y$ and independent variables $x$. The core method is the **method of least squares**: finding a set of parameters that minimizes the sum of squared residuals from all observations to the fitted line.

### 1.2 Formula

$$y_i = \beta_0 + \beta_1 x_{i1} + \cdots + \beta_p x_{ip} + \epsilon_i \qquad \Leftrightarrow \qquad Y = X\beta + \epsilon$$

- **$y_i$**: the value of the dependent variable for the $i$-th observation
- **$\beta_0$**: the intercept term
- **$\beta_1, \ldots, \beta_p$**: partial regression coefficients for the independent variables
- **$x_{i1}, \ldots, x_{ip}$**: values of each independent variable for the $i$-th observation
- **$\epsilon_i$**: the error term (residual)

### 1.3 Core Assumptions (LINE)

The validity of OLS depends on the **Gauss-Markov assumptions**, remembered by the mnemonic **LINE**:

- **L — Linearity**: $E[y|x]$ is strictly linear in $x$
- **I — Independence**: the $\epsilon_i$ are free of autocorrelation
- **N — Normality**: $\epsilon \sim N(0, \sigma^2)$; particularly critical for hypothesis testing (t/F tests)
- **E — Equal variance (Homoscedasticity)**: $\text{Var}(\epsilon_i) = \sigma^2$, constant across all values of $x$
- **(Additional) No multicollinearity**: the independent variables are not highly correlated with one another; severe collinearity inflates parameter variances dramatically

### 1.4 Heteroscedasticity

> When the homoscedasticity assumption is violated, the spread of the errors grows or shrinks as the independent variable changes (the most common pattern is a "fan-shaped" or "funnel-shaped" dispersion).

**Sources:**

- **Large scale differences in the data**: as income rises, the variance of consumption expenditure also rises (lower-income groups have relatively fixed expenditure)
- **Omitted important variables**: the influence of omitted variables is absorbed into the residuals; if those omitted variables are correlated with a predictor, the residual variance changes systematically
- **Misspecified functional form**: when the true relationship is exponential or follows a power law but is forced into a linear fit, the bias increases geometrically with $x$
- **Outliers**: extreme observations pull the regression line and disrupt local variance homogeneity
- **Changing measurement error**: larger objects of measurement produce larger absolute errors (relative error may be stable, but residuals are not)

**Consequences:**

- ✅ **Coefficient estimates remain unbiased**: $\hat{\beta}$ still points in the correct direction; the general trend is preserved
- ❌ **OLS is no longer BLUE**: it loses the minimum-variance property and is no longer the Best Linear Unbiased Estimator
- ❌ **Standard errors are incorrectly estimated**: the conventional OLS formula severely underestimates the SE of the coefficients
- ❌ **Hypothesis tests are completely invalidated**: $t = \hat{\beta}/SE$ is artificially inflated, $p$-values are too small, the risk of **Type I Error (false positives)** increases sharply, and confidence intervals are spuriously narrow

---

## 2. Standard ANOVA

### 2.1 Principles

ANOVA is used to test whether three or more population means are all equal. **In essence, ANOVA is a special case of OLS regression in which all independent variables are categorical.**

The central idea is to partition the total variation into two components and use the **F-statistic** (the ratio of the two) to determine whether the between-group differences are significantly larger than random error:

- **Between-group variance**: variation attributable to different treatments or categories
- **Within-group variance (Error)**: variation attributable to random error

### 2.2 Formula

One-way ANOVA:

$$y_{ij} = \mu + \alpha_i + \epsilon_{ij}$$

- **$y_{ij}$**: the $j$-th observation in group (treatment) $i$
- **$\mu$**: the overall grand mean across all groups
- **$\alpha_i$**: the treatment effect for group $i$ ($= \mu_i - \mu$)
- **$\epsilon_{ij}$**: the random error term

### 2.3 Core Assumptions

- **Independence**: samples are drawn randomly; observations within and across groups are mutually independent
- **Normality**: data within each group follow a normal distribution (equivalently, residuals are normally distributed)
- **Homogeneity of Variance**: the population variances of all comparison groups must be equal (the homoscedasticity assumption of OLS expressed in the context of categorical predictors); if variances are unequal, use **Welch's ANOVA** instead

---

## 3. Gamma Distribution and Gamma GLM

### 3.1 Gamma vs. Poisson: Two Sides of the Same Coin

> **Key distinction**: Poisson describes **discrete integers** (counting events), while Gamma describes **continuous positive real numbers** (measuring time or monetary amounts). Both arise from the same underlying stochastic process — the **Poisson process**.

Using a bus stop as an illustration:

- **Poisson distribution (asking "how many?")**: fix a 1-hour window and count how many buses arrive during that period → the result can only be an integer: 0, 1, 2, 3, …
- **Gamma distribution (asking "how long?")**: start a stopwatch and measure the total time until the $k$-th bus arrives → the result is a continuous positive real number such as 15.2 minutes or 43.78 minutes

**Mathematical relationship**: if events occur according to a Poisson process, the time required to wait for multiple events follows a Gamma distribution. (Special case: the waiting time for a single event follows an **Exponential distribution**, which is a special case of the Gamma.)

### 3.2 When to Use Gamma GLM

Gamma GLM is designed for dependent variables $y$ that simultaneously satisfy all three of the following conditions:

- **Continuous numeric values**: not integer counts, but continuous quantities that can take decimal values
- **Strictly positive**: $y > 0$; values cannot be negative, and typically cannot be zero
- **Long-tailed/right-skewed and heteroscedastic**: the majority of values are concentrated in a small range while a small number of extreme large values produce a long right tail; moreover, **as the mean increases, the variance grows quadratically** (groups with larger means exhibit far greater internal variability)

**Classic application domains:**

- **Insurance claim amounts**: most claims are small (minor scrapes and dents), but occasional catastrophic claims generate an extreme right tail; furthermore, insurance lines with higher average payouts tend to have greater volatility
- **Medical expenditure or length of hospital stay**: most patients incur modest costs and have short stays, while a small number of severely ill patients incur extremely high costs over extended periods

### 3.3 Formula and Link Function

The expected value (mean) of the Gamma distribution is denoted $\mu$. In Gamma GLM:

- **Canonical Link — Inverse (Reciprocal) Link**:

$$\frac{1}{\mu} = \beta_0 + \beta_1 x_1 + \cdots + \beta_p x_p$$

- **The overwhelmingly preferred choice in practice — Log Link** (ensures the predicted mean $\mu$ is always positive; shares the same link function as Poisson regression, which is one reason the two are sometimes confused):

$$\ln(\mu) = \beta_0 + \beta_1 x_1 + \cdots + \beta_p x_p$$

---

## 4. Deep Dive: Link Functions in GLM

### 4.1 Why Do Both Gamma and Poisson Use the Log Link?

The target quantity being predicted (the mean $\mu$) in both Gamma and Poisson models shares a common strict constraint: **it must be greater than zero**. The Log Link is the most natural way to enforce this:

$$\ln(\mu) = X\beta \quad \Rightarrow \quad \mu = e^{X\beta}$$

No matter what positive or negative value $X\beta$ takes, the predicted value $\mu = e^{X\beta}$ **is always greater than zero**.

**Poisson link:**

- **Standard link — Log**: $\ln(\mu) = X\beta$
- **Why?** $X\beta$ can take any real value; exponentiating both sides guarantees $\mu = e^{X\beta} > 0$, perfectly consistent with the non-negativity requirement of count data

**Gamma link:**

- **Theoretical canonical link — Inverse (Reciprocal)**: $1/\mu = X\beta$; rarely used in practice because if $X\beta$ is negative the implied mean becomes negative, violating the Gamma distribution's requirement that values be strictly positive
- **Dominant choice in practice — Log**: $\ln(\mu) = X\beta$; ensures predicted values are always positive; additionally, the Log Link converts additive effects into **multiplicative effects (percentage changes)**, which is substantively meaningful for interpretation (see the blood biomarker discussion below)

### 4.2 Why Does Substituting Poisson for Logistic Regression Yield More Meaningful Risk Ratios?

> This is a classic technique in epidemiology and medical statistics, commonly known as **Modified Poisson Regression**.

#### The "Pain Point" of Logistic Regression: OR Severely Overestimates RR at High Prevalence

Logistic regression inherently produces only OR (Odds Ratio), whereas clinicians and the general public think in terms of **RR (Relative Risk)**:

- When a disease is **very rare (incidence < 10%)**, the numerical difference between OR and RR is negligible
- When a disease is **common (incidence > 10%)**, OR **severely overestimates** the true RR. For example:
  - True RR (Relative Risk) = 80% / 40% = **2.0** (cure rate doubled — easy to understand)
  - Logistic regression gives OR = (0.8/0.2) / (0.4/0.6) = 4 / 0.67 = **6.0**
  - Reporting that "the medication increased the odds of cure by 6-fold" is misleading to readers

#### How Does Poisson Solve This Problem?

By fitting binary outcome data (0 and 1) through a generalized linear model with a Log Link:

$$\ln(p) = \beta_0 + \beta_1 x$$

Exponentiating both sides: $p = e^{\beta_0} \cdot e^{\beta_1 x}$. The exponentiated coefficient $e^{\beta_1}$ for predictor $x$ is **directly interpretable as the RR (Risk Ratio)** in a statistical sense.

> **Note**: Because binary data do not satisfy the variance assumptions of the Poisson distribution, running a straight Poisson regression on 0/1 outcomes will produce p-values and confidence intervals that are too narrow. Therefore, **Robust Standard Errors (sandwich variance estimator)** must be specified in the software. This is the classic solution proposed by Zou (2004).

### 4.3 Why Use Log-Gamma When Studying Blood Biomarkers?

**Blood biomarkers fit the three defining characteristics of a Gamma distribution perfectly:**

- **Strictly non-negative and highly right-skewed**: blood biomarker concentrations in healthy individuals are typically maintained at very low levels (e.g., CRP at 1–3 mg/L), but when inflammation or illness occurs, concentrations do not increase by a few points — they surge tens or hundreds of times over (e.g., spiking to 100–300 mg/L). The data are highly concentrated on the left with an extremely long right tail, which is catastrophic for ordinary OLS
- **Variance increases sharply as the mean increases ($\text{Var} \propto \mu^2$)**: CRP concentrations in healthy individuals cluster around 1–3 mg/L with very little individual variation (small variance); but in a severe infection ward, one patient's CRP may be 50 while another's is 250, with enormous inter-individual variability (large variance). As the mean grows, variance scales quadratically — this is precisely where the Gamma distribution excels
- **Why "Log"-Gamma (multiplicative effects align with biological intuition)**: in a Gamma regression with a Log Link ($\ln(\mu) = X\beta$), the coefficients represent **multiplicative effects (percentage changes)**:
  - Under an ordinary additive model: the interpretation would be "each additional year of age increases biomarker concentration by 2 mg/L." This is biologically implausible — an increase of 2 from a baseline of 1 is a doubling, whereas an increase of 2 from a baseline of 100 is negligible
  - Under a Log-Gamma model: the coefficient is interpreted as "each additional year of age increases biomarker concentration by 5%." In human biochemical metabolism, enzyme kinetics, and pharmacokinetics, concentration changes almost universally follow proportional (multiplicative) patterns

> **Summary**: When studying blood biomarkers, the Gamma distribution precisely captures the long-tailed and heteroscedastic nature of the data, the Log link function guarantees that predicted concentrations are always positive, and the resulting coefficients provide "percentage/fold-change" interpretations that align with human physiology. This is why Log-Gamma is the gold standard for this type of data.

---

## 5. What Does "Linear" Mean in GLM?

**Answer**: In terms of visual appearance in the original data space, GLMs are nonlinear (curved); but at the core of the mathematical equations, they are purely linear models.

- **Visually nonlinear (in the natural data space)**: plotting logistic regression (an S-shaped curve) or Poisson regression (an exponential growth curve) on coordinate axes reveals that the relationship between the dependent variable $Y$ and the independent variable $X$ is unambiguously curved. In terms of the final predictions, these models accommodate nonlinear relationships

- **Mathematically linear (in the parameter space)**: why, then, are they still called "linear" models? Because on the right-hand side of the equation, **the parameters ($\beta$) and the independent variables ($X$) are combined in a strictly additive linear fashion**:

$$\text{Linear predictor} = \beta_0 + \beta_1 X_1 + \beta_2 X_2 + \ldots$$

The term "generalized linear" means that although the model uses a **link function** to "straighten out" the curved real-world data (e.g., by taking the Log or Logit), in the transformed space — after the data have been linearized — the core equation being fitted is still the familiar straight-line equation.

> **Summary**: A GLM is a "linear" model in "nonlinear" clothing. It cleverly uses a mathematical transformation (the link function) to retain the computational simplicity and interpretive clarity of linear equations, while simultaneously solving the problem of fitting nonlinear data.

---

## 6. Is Modified Poisson Regression Necessary for Rare Binary Outcomes?

**Conclusion**: it is technically feasible in software, but from both a statistical and practical standpoint — **it is entirely unnecessary and is in fact not recommended**.

### 6.1 The Rare Disease Assumption

A fundamental principle in epidemiology states: **when the event rate is low, the OR converges mathematically to the RR ($OR \approx RR$)**.

- **Derivation**: $\text{Odds} = \frac{p}{1-p}$. If the incidence $p$ is very small (e.g., 0.01), then the denominator $1-p$ is extremely close to 1, so $\text{Odds} \approx p$. Consequently, the ratio of two odds (OR) approximates the ratio of two probabilities (RR)
- **Practical implication**: since OR equals RR when outcomes are rare, one can simply apply the classic, well-established Logistic regression to obtain an OR and interpret it directly as an RR when communicating with clinicians — no need to take the roundabout route of Modified Poisson

### 6.2 Potential Pitfalls of Modified Poisson

Forcing Modified Poisson on rare outcomes may introduce problems:

- **Convergence issues**: Logistic regression is purpose-built for binary data; its maximum likelihood estimation (MLE) algorithm is highly stable when processing sparse 0/1 data. Poisson regression is fundamentally designed for count data (0, 1, 2, 3, …); forcing it to handle data consisting almost entirely of zeros (with very few ones) may cause certain statistical software packages to produce errors or fail to converge
- **Boundary overflow risk**: the Logit link in Logistic regression acts as a built-in "safety lock," guaranteeing that all predicted probabilities remain within [0, 1]. The Log link in Poisson regression has no upper bound; although predicted values are unlikely to exceed 1 for rare events, in extreme scenarios (or when some predictors take extreme values), it is mathematically possible for the model to predict an "event probability > 1 (i.e., > 100%)," which is an absurd result

### 6.3 Practical Recommendations

Before analyzing binary outcome data, first examine the proportion of events (Y = 1):

- **Event rate > 10%**: use **Modified Poisson** (with Robust standard errors) to estimate RR
- **Event rate < 10%**: use the straightforward **Logistic regression**, compute the OR, and interpret it confidently as an approximation of RR

---

## 7. OLS Closed-Form Solution vs. GLM Iterative Numerical Solution

### 7.1 OLS: A One-Step Closed-Form Solution

The sole objective of OLS is to minimize the sum of squared vertical distances (residuals) from all observations to the regression line.

Because OLS assumes a continuous dependent variable with constant variance (homoscedasticity), its mathematical formulation is elegant. In matrix algebra, minimizing the sum of squared residuals is equivalent to finding the minimum of a paraboloid — one simply differentiates with respect to the parameters, sets the derivative to zero, and solves. The exact answer can be **computed directly from a single fixed formula, with no iteration required**:

$$\hat{\beta} = (X^T X)^{-1} X^T Y$$

> **Note**: as long as there is no perfect multicollinearity among the predictors (i.e., the matrix $X^TX$ is invertible), the computer performs a single matrix multiplication and inversion to obtain the exact value of $\beta$ instantaneously.

### 7.2 GLM: Maximum Likelihood Estimation and Iterative Methods (MLE + IRLS)

Generalized linear models (such as Logistic, Poisson, and Gamma) deal with data that are non-normally distributed, involve nonlinear link functions, and have variances that change with the mean. In this setting, the elegant closed-form formula of least squares breaks down.

#### Core Idea — Maximum Likelihood Estimation (MLE)

Rather than minimizing the sum of squared errors, GLM seeks to **maximize the probability (likelihood) of observing the data at hand**. It asks: for what value of the parameter $\beta$ would the observed data be most likely to occur?

#### Solution Method — Iteratively Reweighted Least Squares (IRLS)

Because complex link functions (such as Logit or Log) are involved, the likelihood function becomes highly complex, and **there is no closed-form formula analogous to OLS**. The computer must employ a numerical optimization algorithm — typically Iteratively Reweighted Least Squares (IRLS), or equivalently the Newton-Raphson method — to search iteratively for the optimum:

- **Step 1 (initial guess)**: the computer assigns an arbitrary starting value to $\beta$ (e.g., all zeros)
- **Step 2 (compute weights and residuals)**: based on the current $\beta$ and the chosen link function, the predicted value for each observation is computed. Because the variance in a GLM changes with the mean, observations with smaller variance (greater certainty) are assigned **larger weights**
- **Step 3 (update direction)**: using the weighted residuals, the gradient of the likelihood function at the current point is computed, and $\beta$ is updated by stepping in the direction that increases the likelihood
- **Step 4 (repeat)**: the new $\beta$ is used to recompute the weights (hence "iteratively reweighted"), and the direction is updated again
- **Step 5 (convergence)**: when successive estimates of $\beta$ differ by a negligible amount (e.g., less than 0.0001), the computer concludes that it has reached the maximum of the likelihood function, declares **model convergence**, halts the computation, and reports the results

---

## 8. Design Matrix

### 8.1 Core Concept

Multiple regression is expressed compactly in matrix form as: $Y' = X\hat{\beta}$

- **$Y'$**: a column vector of predicted values
- **$X$ (Design Matrix)**: a matrix containing a column of ones (for the intercept) and columns of observed values for all predictors
- **$\hat{\beta}$**: a vector of estimated regression coefficients

**Expanded form with a single predictor (5 observations):**

$$Y'_i = \hat{\beta}_0 + \hat{\beta}_1 X_{1i} \;\equiv\; \hat{\beta}_0(1) + \hat{\beta}_1 X_{1i}$$

$$\begin{bmatrix} Y'_1 \\ Y'_2 \\ Y'_3 \\ Y'_4 \\ Y'_5 \end{bmatrix} = \begin{bmatrix} 1 & X_{11} \\ 1 & X_{12} \\ 1 & X_{13} \\ 1 & X_{14} \\ 1 & X_{15} \end{bmatrix} \begin{bmatrix} \hat{\beta}_0 \\ \hat{\beta}_1 \end{bmatrix} = \begin{bmatrix} \hat{\beta}_0 + \hat{\beta}_1 X_{11} \\ \hat{\beta}_0 + \hat{\beta}_1 X_{12} \\ \hat{\beta}_0 + \hat{\beta}_1 X_{13} \\ \hat{\beta}_0 + \hat{\beta}_1 X_{14} \\ \hat{\beta}_0 + \hat{\beta}_1 X_{15} \end{bmatrix}$$

**Expanded form with two predictors:**

$$Y'_i = \hat{\beta}_0 + \hat{\beta}_1 X_{1i} + \hat{\beta}_2 X_{2i}$$

$$\begin{bmatrix} Y'_1 \\ Y'_2 \\ Y'_3 \\ Y'_4 \\ Y'_5 \end{bmatrix} = \begin{bmatrix} 1 & X_{11} & X_{21} \\ 1 & X_{12} & X_{22} \\ 1 & X_{13} & X_{23} \\ 1 & X_{14} & X_{24} \\ 1 & X_{15} & X_{25} \end{bmatrix} \begin{bmatrix} \hat{\beta}_0 \\ \hat{\beta}_1 \\ \hat{\beta}_2 \end{bmatrix} = \begin{bmatrix} \hat{\beta}_0 + \hat{\beta}_1 X_{11} + \hat{\beta}_2 X_{21} \\ \hat{\beta}_0 + \hat{\beta}_1 X_{12} + \hat{\beta}_2 X_{22} \\ \hat{\beta}_0 + \hat{\beta}_1 X_{13} + \hat{\beta}_2 X_{23} \\ \hat{\beta}_0 + \hat{\beta}_1 X_{14} + \hat{\beta}_2 X_{24} \\ \hat{\beta}_0 + \hat{\beta}_1 X_{15} + \hat{\beta}_2 X_{25} \end{bmatrix}$$

### 8.2 OLS Matrix Solution

The coefficient vector is obtained in one step via matrix algebra: construct the design matrix $X$, compute its transpose $X'$, form the product $X'X$, invert it, and multiply by $X'Y$:

$$\hat{\beta} = (X'X)^{-1}X'Y$$

> As long as $X'X$ is invertible (no perfect multicollinearity), this formula yields the exact solution directly, with no iteration.

### 8.3 Covariance Matrix of the Coefficient Estimates

The covariance matrix of the vector $\hat{\beta}$:

$$\text{cov}(\hat{\beta}) = \sigma^2 (X'X)^{-1}$$

where the residual variance $\sigma^2$ is estimated from the residual vector $\epsilon$:

$$\sigma^2 = \frac{\epsilon'\epsilon}{n - s}$$

- **$\epsilon$**: the residual vector ($= Y - X\hat{\beta}$)
- **$s$**: the number of regression coefficients estimated in the model (including the intercept)
- **$n - s$**: the residual degrees of freedom

The diagonal elements of the covariance matrix are the variances of the individual coefficient estimates; their square roots are the **standard errors (SE)** of the respective coefficients. This matrix can be extracted in software via the `vcov()` function to verify the SE values reported in the output.

### 8.4 Matrix Representation of $\text{Var}(\mathbf{y})$ in OLS

The OLS data-generating model is $\mathbf{y} = X\beta + \boldsymbol{\epsilon}$, where $X\beta$ is fixed (non-random), so:

$$\text{Var}(\mathbf{y}) = \text{Var}(\boldsymbol{\epsilon}) = \sigma^2 \mathbf{I}_n$$

Expanded as an $n \times n$ matrix:

$$\text{Var}(\mathbf{y}) = \begin{bmatrix} \sigma^2 & 0 & \cdots & 0 \\ 0 & \sigma^2 & \cdots & 0 \\ \vdots & & \ddots & \vdots \\ 0 & 0 & \cdots & \sigma^2 \end{bmatrix}$$

**The diagonal elements are the variance of each observation (all equal); the off-diagonal elements are all zero, indicating mutual independence across observations.**

This is precisely the two core assumptions of OLS expressed in the language of covariance matrices: **Homoscedasticity + Independence**.

#### Correspondence to MLM / LSEM

- **OLS**: $\text{Var}(\mathbf{y}) = \sigma^2 \mathbf{I}_n$ — no random effects; residuals are independent and identically distributed
- **MLM**: $V_i = Z_i G Z_i' + R_i$ — random effects plus within-individual residual structure
- **LSEM**: $\Sigma(\theta) = \Lambda\Psi\Lambda' + \Theta$ — latent factors plus indicator residuals

OLS is the **most constrained special case** of the three: setting $Z_i = \mathbf{0}$ (no random effects) and $R_i = \sigma^2 \mathbf{I}$ (homogeneous independent residuals) reduces the MLM formula immediately to the OLS form $\sigma^2\mathbf{I}$; likewise, setting $\Lambda = \mathbf{0}$ and $\Theta = \sigma^2\mathbf{I}$ reduces the LSEM to OLS.

#### Extensions When Assumptions Are Violated

- **Heteroscedasticity**: off-diagonal elements remain zero (independence holds), but diagonal elements differ:

$$\text{Var}(\mathbf{y}) = \begin{bmatrix} \sigma_1^2 & & \\ & \ddots & \\ & & \sigma_n^2 \end{bmatrix} = \Omega$$

- **Autocorrelation (e.g., AR(1))**: diagonal elements are equal, but off-diagonal elements are non-zero:

$$\text{Var}(\mathbf{y}) = \sigma^2 \begin{bmatrix} 1 & \rho & \rho^2 & \cdots \\ \rho & 1 & \rho & \cdots \\ \rho^2 & \rho & 1 & \cdots \\ \vdots & & & \ddots \end{bmatrix}$$

- **General case (GLS)**: $\text{Var}(\mathbf{y}) = \sigma^2\Omega$; the OLS formula becomes the Generalized Least Squares (GLS) estimator:

$$\hat{\beta}_{\text{GLS}} = (X'\Omega^{-1}X)^{-1}X'\Omega^{-1}Y$$

When $\Omega = \mathbf{I}$, this reduces back to the OLS estimator $\hat{\beta} = (X'X)^{-1}X'Y$.

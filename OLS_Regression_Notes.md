# Statistical Modeling Notes

---

## Ordinary Least Squares Regression (OLS Regression)

### Principles

OLS regression examines the linear relationship between a continuous dependent variable $y$ and independent variables $x$. The core method is the **method of least squares**: finding a set of parameters that minimizes the sum of squared residuals from all observations to the fitted line.

### Formula

$$
y_i = \beta_0 + \beta_1 x_{i1} + \cdots + \beta_p x_{ip} + \epsilon_i \qquad \Leftrightarrow \qquad Y = X\beta + \epsilon
$$

- **$y_i$**: the value of the dependent variable for the $i$-th observation
- **$\beta_0$**: the intercept term
- **$\beta_1, \ldots, \beta_p$**: partial regression coefficients for the independent variables
- **$x_{i1}, \ldots, x_{ip}$**: values of each independent variable for the $i$-th observation
- **$\epsilon_i$**: the error term (residual)

### Core Assumptions (LINE)

The validity of OLS depends on the **Gauss-Markov assumptions**, remembered by the mnemonic **LINE**:

- **L — Linearity**: $E[y|x]$ is strictly linear in $x$
- **I — Independence**: the $\epsilon_i$ are free of autocorrelation
- **N — Normality**: $\epsilon \sim N(0, \sigma^2)$; particularly critical for hypothesis testing (t/F tests)
- **E — Equal variance (Homoscedasticity)**: $\text{Var}(\epsilon_i) = \sigma^2$, constant across all values of $x$
- **(Additional) No multicollinearity**: the independent variables are not highly correlated with one another; severe collinearity inflates parameter variances dramatically

### Heteroscedasticity

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

## Standard ANOVA

### Principles

ANOVA is used to test whether three or more population means are all equal. **In essence, ANOVA is a special case of OLS regression in which all independent variables are categorical.**

The central idea is to partition the total variation into two components and use the **F-statistic** (the ratio of the two) to determine whether the between-group differences are significantly larger than random error:

- **Between-group variance**: variation attributable to different treatments or categories
- **Within-group variance (Error)**: variation attributable to random error

### Formula

One-way ANOVA:

$$
y_{ij} = \mu + \alpha_i + \epsilon_{ij}
$$

- **$y_{ij}$**: the $j$-th observation in group (treatment) $i$
- **$\mu$**: the overall grand mean across all groups
- **$\alpha_i$**: the treatment effect for group $i$ ($= \mu_i - \mu$)
- **$\epsilon_{ij}$**: the random error term

### Core Assumptions

- **Independence**: samples are drawn randomly; observations within and across groups are mutually independent
- **Normality**: data within each group follow a normal distribution (equivalently, residuals are normally distributed)
- **Homogeneity of Variance**: the population variances of all comparison groups must be equal (the homoscedasticity assumption of OLS expressed in the context of categorical predictors); if variances are unequal, use **Welch's ANOVA** instead

---

## What Does "Linear" Mean in GLM?

**Answer**: In terms of visual appearance in the original data space, GLMs are nonlinear (curved); but at the core of the mathematical equations, they are purely linear models.

- **Visually nonlinear (in the natural data space)**: plotting logistic regression (an S-shaped curve) or Poisson regression (an exponential growth curve) on coordinate axes reveals that the relationship between the dependent variable $Y$ and the independent variable $X$ is unambiguously curved. In terms of the final predictions, these models accommodate nonlinear relationships
- **Mathematically linear (in the parameter space)**: why, then, are they still called "linear" models? Because on the right-hand side of the equation, **the parameters ($\beta$) and the independent variables ($X$) are combined in a strictly additive linear fashion**:

$$
\text{Linear predictor} = \beta_0 + \beta_1 X_1 + \beta_2 X_2 + \ldots
$$

The term "generalized linear" means that although the model uses a **link function** to "straighten out" the curved real-world data (e.g., by taking the Log or Logit), in the transformed space — after the data have been linearized — the core equation being fitted is still the familiar straight-line equation.

> **Summary**: A GLM is a "linear" model in "nonlinear" clothing. It cleverly uses a mathematical transformation (the link function) to retain the computational simplicity and interpretive clarity of linear equations, while simultaneously solving the problem of fitting nonlinear data.

---

## OLS Closed-Form Solution vs. GLM Iterative Numerical Solution

### OLS: A One-Step Closed-Form Solution

The sole objective of OLS is to minimize the sum of squared vertical distances (residuals) from all observations to the regression line.

Because OLS assumes a continuous dependent variable with constant variance (homoscedasticity), its mathematical formulation is elegant. In matrix algebra, minimizing the sum of squared residuals is equivalent to finding the minimum of a paraboloid — one simply differentiates with respect to the parameters, sets the derivative to zero, and solves. The exact answer can be **computed directly from a single fixed formula, with no iteration required**:

$$
\hat{\beta} = (X^T X)^{-1} X^T Y
$$

> **Note**: as long as there is no perfect multicollinearity among the predictors (i.e., the matrix $X^TX$ is invertible), the computer performs a single matrix multiplication and inversion to obtain the exact value of $\beta$ instantaneously.

### GLM: Maximum Likelihood Estimation and Iterative Methods (MLE + IRLS)

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

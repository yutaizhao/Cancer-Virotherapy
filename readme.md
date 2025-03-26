Up to date model (Bajzer) :

$$
\begin{aligned}
y^{\prime} & =r y\left[1-(y+x)^{\epsilon}/ K^{\epsilon}\right]-\kappa y v-\rho x y \\
x^{\prime} & =\kappa y v-\delta x \\
v^{\prime} & =\alpha x-\omega v-\kappa y v
\end{aligned}
$$
  
We improved it :

In this model, we separate the syncytia population (s) from the infected cell population (x). Because syncytia should not be able to release virus nor infect uninfected cells.

We also take into account the fact that once intertumoral infection occurs (represented by $\rho x y$), uninfected cells (y) and infected cells (x) interact to form syncytia. As a result, the populations of both infected and uninfected cells should decrease, while the syncytia population should increase by the same amount. In other words, the term $+\rho x y$ should not be neglected. 

Finally, we introduce the parameters $\gamma$ death rate of syncytia and $l \gamma s$ to represent the potential viral release when syncytia die.
 
$$
\begin{aligned}
y^{\prime} & =r y\left[1-(y+x)^{\epsilon} / K^{\epsilon}\right]-\kappa y v-\rho x y \\
x^{\prime} & =\kappa y v-\delta x - \rho x y \\
s^{\prime} & = 2 \rho x y - \gamma s \\
v^{\prime} & =\alpha x-\omega v-\kappa y v + l \gamma s
\end{aligned}
 $$  

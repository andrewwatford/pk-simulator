### Equations
---


$$
\frac{d q_{0}}{d t} = D_{0}(t)- k_{0,1}\frac{q_{0}}{V_{0}}
$$



$$
\frac{d q_{1}}{d t} = - C_{1}\frac{q_{1}}{V_{1}}+ k_{0,1}\frac{q_{0}}{V_{0}}- k_{1,2}\left(\frac{q_{1}}{V_{1}} - \frac{q_{2}}{V_{2}}\right)
$$



$$
\frac{d q_{2}}{d t} = - k_{1,2}\left(\frac{q_{2}}{V_{2}} - \frac{q_{1}}{V_{1}}\right)
$$

### Compartments
---


| Index | Compartment |
|-------|-------------|
| 0     | absorbing     |
| 1     | central     |
| 2     | peripheral     |
### Variable definitions
---


| Symbol | Quantity |
|-------|-------------|
| $t$ | Time |
| $q_i$ | Mass of drug in compartment i |
| $V_i$ | Volume of compartment i |
| $D_i$ | Dosage into compartment i |
| $C_i$ | Clearance rate from compartment i |
| $k_{i,j}$ | Rate constant for flux between compartments i and j |

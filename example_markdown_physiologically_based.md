### Equations
---


$$
\frac{d q_{0}}{d t} = D_{0}(t)- k_{0,1}\frac{q_{0}}{V_{0}}
$$



$$
\frac{d q_{1}}{d t} = + k_{0,1}\frac{q_{0}}{V_{0}}- k_{3,1}\left(\frac{q_{1}}{V_{1}} - \frac{q_{3}}{V_{3}}\right)- k_{2,1}\left(\frac{q_{1}}{V_{1}} - \frac{q_{2}}{V_{2}}\right)
$$



$$
\frac{d q_{2}}{d t} = - k_{2,1}\left(\frac{q_{2}}{V_{2}} - \frac{q_{1}}{V_{1}}\right)- k_{2,3}\left(\frac{q_{2}}{V_{2}} - \frac{q_{3}}{V_{3}}\right)- k_{2,4}\left(\frac{q_{2}}{V_{2}} - \frac{q_{4}}{V_{4}}\right)- k_{2,5}\left(\frac{q_{2}}{V_{2}} - \frac{q_{5}}{V_{5}}\right)- k_{2,6}\left(\frac{q_{2}}{V_{2}} - \frac{q_{6}}{V_{6}}\right)- k_{2,7}\left(\frac{q_{2}}{V_{2}} - \frac{q_{7}}{V_{7}}\right)- k_{2,8}\left(\frac{q_{2}}{V_{2}} - \frac{q_{8}}{V_{8}}\right)
$$



$$
\frac{d q_{3}}{d t} = - k_{3,1}\left(\frac{q_{3}}{V_{3}} - \frac{q_{1}}{V_{1}}\right)- k_{3,4}\left(\frac{q_{3}}{V_{3}} - \frac{q_{4}}{V_{4}}\right)- k_{2,3}\left(\frac{q_{3}}{V_{3}} - \frac{q_{2}}{V_{2}}\right)
$$



$$
\frac{d q_{4}}{d t} = - k_{3,4}\left(\frac{q_{4}}{V_{4}} - \frac{q_{3}}{V_{3}}\right)- k_{2,4}\left(\frac{q_{4}}{V_{4}} - \frac{q_{2}}{V_{2}}\right)
$$



$$
\frac{d q_{5}}{d t} = - k_{2,5}\left(\frac{q_{5}}{V_{5}} - \frac{q_{2}}{V_{2}}\right)
$$



$$
\frac{d q_{6}}{d t} = - k_{2,6}\left(\frac{q_{6}}{V_{6}} - \frac{q_{2}}{V_{2}}\right)
$$



$$
\frac{d q_{7}}{d t} = - k_{2,7}\left(\frac{q_{7}}{V_{7}} - \frac{q_{2}}{V_{2}}\right)
$$



$$
\frac{d q_{8}}{d t} = - k_{2,8}\left(\frac{q_{8}}{V_{8}} - \frac{q_{2}}{V_{2}}\right)
$$

### Compartments
---


| Index | Compartment |
|-------|-------------|
| 0     | ingestion     |
| 1     | gut     |
| 2     | arteries     |
| 3     | liver     |
| 4     | pancreas     |
| 5     | heart_and_lung     |
| 6     | brain     |
| 7     | muscle_and_skin     |
| 8     | kidney     |
### Variable definitions
---


| Symbol | Quantity |
|-------|-------------|
| $t$ | Time |
| $q_i$ | Mass of drug in compartment i |
| $V_i$ | Volume of compartment i |
| $D_i$ | Dosage into compartment i |
| $k_{i,j}$ | Rate constant for flux between compartments i and j |

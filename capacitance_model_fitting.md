## Capacitance vs reverse bias

$$C_{j} = \frac{A}{2} \sqrt{\frac{2 q \varepsilon N_{D}^{+}(T)}{V_0 (T) - V}} + C_{p}.$$

Truncate terms, otherwise scipy.curve_fit() gives non-physical values.
$$a = \frac{A}{2} \sqrt{2 q \varepsilon N_{D}^{+}(T)},$$
$$C_{j, \text{trunc}} = \frac{a}{\sqrt{V_0 (T) - V}} + C_{p}.$$

$C_{p}$ accounts for parasitic capacitance in the experimental set up.

## Capacitance vs temperature
$$C(T) = \frac{A}{2} \sqrt{\frac{2 q \varepsilon N_{D}^{+}(T)}{V_0 (T) - V}} + C_{p}.$$

Define variable $a,$
$$a = \frac{A}{2} \sqrt{2 q \varepsilon},$$
$$C(T) = a \sqrt{\frac{N_{D}^{+}(T)}{V_{0}(T)-V}} + C_{p},$$

Built-in voltage as a function of T (and ionized acceptors/donors concentration),
$$V_{0}(T) = \frac{k T}{q} \log\left[\frac{N_{D}^{+}(T) N_{A}^{-}(T)}{4} \frac{h^{6}}{(2 \pi k T)^3} (m_{e} m_{h})^{-3/2}\right] + \frac{E_{G}}{q},$$

Expand log of products to sum of logs,
$$V_{0}(T) = \frac{k T}{q} \log\left[\frac{h^{6}}{4 (2 \pi k)^3} (m_{e} m_{h})^{-3/2}\right] - \frac{3 k T}{q} \log T + \frac{k T}{q} \log[N_{D}^{+}(T) N_{A}^{-}(T)] + \frac{E_{G}}{q},$$

Ionized dopants as a function T,
$$N_{D}^{+}(T) = \frac{N_{D}}{\exp\left[\frac{E_{\text{ion, D}}}{2 k T}\right] + 1},$$
$$\log N_{D}^{+}(T) = \log N_{D} - \log\left(\exp\left[\frac{E_{\text{ion, D}}}{2 k T}\right] + 1\right) \approx \log N_{D} - \frac{E_{\text{ion, D}}}{2 k T},$$
$$N_{A}^{-}(T) = \frac{N_{A}}{\exp\left[\frac{E_{\text{ion, A}}}{2 k T}\right] + 1},$$
$$\log N_{A}^{-}(T) = \log N_{A} - \log\left(\exp\left[\frac{E_{\text{ion, A}}}{2 k T}\right] + 1\right) \approx \log N_{A} - \frac{E_{\text{ion, A}}}{2 k T},$$
Then,
$$\log[N_{D}^{+}(T) N_{A}^{-}(T)] = \log N_{D}^{+}(T) + \log N_{A}^{-}(T) \approx \log N_{D} + \log N_{A} - \left[\frac{E_{\text{ion, D}}}{2 k T} + \frac{E_{\text{ion, A}}}{2 k T}\right]$$ 

Plug the above expression back to $V_{0}(T),$
$$V_{0}(T) \approx \frac{k T}{q} \log\left[\frac{h^{6}}{4 (2 \pi k)^3} (m_{e} m_{h})^{-3/2}\right] - \frac{3 k T}{q} \log T + \frac{k T}{q} \log(N_{D}N_{A}) - \frac{E_{\text{ion, D}} + E_{\text{ion, A}}}{2 q} + \frac{E_{G}}{q},$$

Reverse bias is fixed, so
$$V_{0}(T) - V \approx \frac{k T}{q} \log\left[\frac{h^{6}}{4 (2 \pi k)^3} (m_{e} m_{h})^{-3/2} N_{D}N_{A}\right] - \frac{3 k T}{q} \log T - \frac{E_{\text{ion, D}} + E_{\text{ion, A}}}{2 q} + \frac{E_{G}}{q} - V,$$

Now, functionally, $V_{0}(T)$ takes the form of $c_{1} T + c_{2} T \log T + c_{3}.$ Plug this expression back to $C_{j}$ yields
$$C(T) = a \sqrt{\frac{N_{D}^{+}(T)}{c_{1} T + c_{2} T \log T + c_{3}}} + C_{p},$$
$$C(T) = a \sqrt{\frac{1}{c_{1} T + c_{2} T \log T + c_{3}} \frac{N_{D}}{\exp\left[\frac{E_{\text{ion, D}}}{2 k T}\right] + 1}} + C_{p},$$

Let $b = \exp\left[\frac{E_{\text{ion, D}}}{2 k}\right],$
$$C(T) = a \sqrt{\frac{1}{c_{1} T + c_{2} T \log T + c_{3}} \frac{N_{D}}{b^{1/T} + 1}} + C_{p},$$

Absorb $N_{D}$ into variable $a$ by re-defining it as $a = \frac{A}{2} \sqrt{2 q \varepsilon N_{D}}.$ Thus,
$$C(T) = a \left[(c_{1} T + c_{2} T \log T + c_{3}) (1 + b^{1/T})\right]^{-1/2} + C_{p}.$$
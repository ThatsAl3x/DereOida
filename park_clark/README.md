# Clarke- und Park-Transformation für Umrichter / MMC-HVDC

Lernskript mit begleitendem Python-Code (`park_clark_demo.py`).
Starten mit `python3 park_clark/park_clark_demo.py`. Das Skript rechnet alle
Formeln hier nach und speichert Plots in `park_clark/plots/`.

---

## 0. Warum das Ganze?

Ein netzgekoppelter Umrichter sieht drei sinusförmige Größen (Spannung,
Strom) mit 50 Hz. Ein PI-Regler kann einem Sinus **nicht** ohne
Regelabweichung folgen (er hat nur bei f = 0 unendliche Verstärkung).

Idee:
1. **Clarke**: 3 Phasen (abc) → 2 orthogonale Achsen (αβ) + Nullsystem (0).
   Immer noch Sinusgrößen, aber nur noch 2 statt 3.
2. **Park**: αβ in ein Koordinatensystem drehen, das **mit der Netzfrequenz
   mitrotiert** (dq). Dann werden aus den 50-Hz-Sinusgrößen
   **Gleichgrößen** → PI-Regler funktionieren perfekt, und
   Wirk-/Blindleistung entkoppeln sich in d und q.

```
 abc ──Clarke──► αβ0 ──Park(θ)──► dq0      (θ kommt von der PLL)
 abc ◄─inv.Clarke─ αβ0 ◄─inv.Park(θ)── dq0  (Sollspannung → Modulator)
```

---

## 1. Raumzeiger (Space Vector) – die Intuition

Fasse die drei Phasengrößen zu einer komplexen Zahl zusammen:

$$
\underline{x} = \tfrac{2}{3}\left(x_a + \underline{a}\,x_b + \underline{a}^2 x_c\right),
\qquad \underline{a} = e^{j\,2\pi/3}
$$

Für ein symmetrisches System $x_a = \hat X\cos(\omega t)$,
$x_b = \hat X\cos(\omega t - 120°)$, $x_c = \hat X\cos(\omega t + 120°)$ gilt:

$$
\underline{x} = \hat X\, e^{j\omega t}
$$

→ ein Zeiger mit **konstanter Länge** $\hat X$, der mit $\omega$ im Kreis
dreht. Clarke = Real- und Imaginärteil dieses Zeigers. Park = denselben
Zeiger aus einem mitdrehenden Koordinatensystem anschauen → er steht still.

---

## 2. Clarke-Transformation (abc → αβ0)

### Amplitudeninvariant (Standard in der Regelungstechnik, Faktor 2/3)

$$
\begin{bmatrix} x_\alpha \\ x_\beta \\ x_0 \end{bmatrix}
= \frac{2}{3}
\begin{bmatrix}
1 & -\tfrac12 & -\tfrac12 \\
0 & \tfrac{\sqrt3}{2} & -\tfrac{\sqrt3}{2} \\
\tfrac12 & \tfrac12 & \tfrac12
\end{bmatrix}
\begin{bmatrix} x_a \\ x_b \\ x_c \end{bmatrix}
$$

Ausgeschrieben:
- $x_\alpha = \tfrac23\left(x_a - \tfrac12 x_b - \tfrac12 x_c\right)$ (bei $x_0=0$: $x_\alpha = x_a$)
- $x_\beta = \tfrac{1}{\sqrt3}(x_b - x_c)$
- $x_0 = \tfrac13(x_a + x_b + x_c)$

Eigenschaft: Amplitude von $x_\alpha, x_\beta$ = Phasenamplitude $\hat X$.
Preis: Leistung braucht den Faktor 3/2 (siehe Abschnitt 4).

### Leistungsinvariant (Faktor √(2/3))

Gleiche Matrix, aber Vorfaktor $\sqrt{2/3}$ und Nullzeile $1/\sqrt2$.
Die Matrix ist dann orthonormal ($T^{-1} = T^\top$) und
$p = v_\alpha i_\alpha + v_\beta i_\beta + v_0 i_0$ ohne 3/2.
Amplituden sind aber um $\sqrt{3/2}$ größer.

> **Merke:** Immer prüfen, welche Variante ein Paper/Simulink-Block
> verwendet – das ist die häufigste Fehlerquelle (Faktor 3/2 bzw. √(3/2)).

### Inverse (amplitudeninvariant)

- $x_a = x_\alpha + x_0$
- $x_b = -\tfrac12 x_\alpha + \tfrac{\sqrt3}{2} x_\beta + x_0$
- $x_c = -\tfrac12 x_\alpha - \tfrac{\sqrt3}{2} x_\beta + x_0$

---

## 3. Park-Transformation (αβ → dq)

Drehung um den Winkel $\theta$ (Winkel der d-Achse gegenüber α):

$$
\begin{bmatrix} x_d \\ x_q \end{bmatrix}
=
\begin{bmatrix}
\cos\theta & \sin\theta \\
-\sin\theta & \cos\theta
\end{bmatrix}
\begin{bmatrix} x_\alpha \\ x_\beta \end{bmatrix}
\qquad\Longleftrightarrow\qquad
\underline{x}_{dq} = \underline{x}_{\alpha\beta}\, e^{-j\theta}
$$

Inverse: $\underline{x}_{\alpha\beta} = \underline{x}_{dq}\,e^{+j\theta}$, also
$x_\alpha = x_d\cos\theta - x_q\sin\theta$, $x_\beta = x_d\sin\theta + x_q\cos\theta$.
Das Nullsystem bleibt unverändert ($x_0$ geht 1:1 durch).

**Beispiel:** $v_a = \hat V\cos(\omega t + \varphi)$, $\theta = \omega t$
→ $v_d = \hat V\cos\varphi$, $v_q = \hat V\sin\varphi$ – **Gleichgrößen**.
Mit einer PLL, die $\theta$ auf den Netzspannungszeiger einrastet, ist
$\varphi = 0$, also $v_d = \hat V$, $v_q = 0$.

### Konventionen (Achtung!)

| Variante | Ausrichtung | Ergebnis für $v_a=\hat V\cos\omega t$ |
|---|---|---|
| cos-basiert, d auf Phase a (hier verwendet) | d-Achse = α bei θ=0 | $v_d=\hat V,\ v_q=0$ |
| sin-basiert / Simulink „90° behind phase A“ | q-Achse auf Phase a | $v_d=0,\ v_q=\hat V$ |

Manche Lehrbücher (Maschinen) legen q **vor** d (wie hier), andere
dahinter → Vorzeichen von $v_q$ und $Q$ drehen sich. Vor jeder Rechnung
die Konvention festhalten.

---

## 4. Leistung in dq

Komplexe Leistung (amplitudeninvariant): $\underline S = \tfrac32\,\underline v\,\underline i^*$

$$
P = \tfrac32\,(v_d i_d + v_q i_q), \qquad
Q = \tfrac32\,(v_q i_d - v_d i_q)
$$

Mit PLL-Ausrichtung ($v_q = 0$):

$$
P = \tfrac32 v_d\, i_d, \qquad Q = -\tfrac32 v_d\, i_q
$$

→ **$i_d$ stellt die Wirkleistung, $i_q$ die Blindleistung** – völlig
entkoppelt. Das ist der Kern der Vektorregelung (Vector Current Control).

---

## 5. Das Streckenmodell in dq – woher die Kopplungsterme kommen

Umrichter (innere Spannung $v_c$) über $R, L$ am Netz ($v_g$):

$$
L\frac{d\underline i_{\alpha\beta}}{dt} = \underline v_{c,\alpha\beta} - \underline v_{g,\alpha\beta} - R\,\underline i_{\alpha\beta}
$$

Einsetzen von $\underline x_{\alpha\beta} = \underline x_{dq}e^{j\theta}$ und
Produktregel: $\frac{d}{dt}(\underline i_{dq} e^{j\theta}) = (\dot{\underline i}_{dq} + j\omega\,\underline i_{dq})e^{j\theta}$:

$$
L\frac{d\underline i_{dq}}{dt} = \underline v_{c} - \underline v_{g} - R\,\underline i_{dq} - j\omega L\,\underline i_{dq}
$$

In Komponenten:

$$
\begin{aligned}
L\frac{di_d}{dt} &= v_{cd} - v_{gd} - R i_d + \omega L\, i_q \\
L\frac{di_q}{dt} &= v_{cq} - v_{gq} - R i_q - \omega L\, i_d
\end{aligned}
$$

Die Terme $\pm\omega L i$ koppeln d und q. Der **Stromregler** kompensiert
sie (Entkopplung + Netzspannungsvorsteuerung):

$$
\begin{aligned}
v_{cd}^* &= v_{gd} + \mathrm{PI}(i_d^* - i_d) - \omega L\, i_q \\
v_{cq}^* &= v_{gq} + \mathrm{PI}(i_q^* - i_q) + \omega L\, i_d
\end{aligned}
$$

Danach sieht jeder Kanal nur noch $L\,di/dt + Ri = \mathrm{PI}(\dots)$ –
eine einfache PT1-Strecke. Typische Auslegung (Internal Model Control):
$K_p = \alpha_c L$, $K_i = \alpha_c R$ mit Bandbreite $\alpha_c$
(z. B. $2\pi\cdot 300$ Hz) → geschlossener Kreis = PT1 mit $\alpha_c$.

### Regelkaskade eines VSC/MMC

```
 P* / Vdc* ──PI──► id* ─┐                         ┌──► v_cd* ─┐
                        ├─► Stromregler (dq) ─────┤            ├─► inv. Park(θ) ► inv. Clarke ► Modulation
 Q* / Vac* ──PI──► iq* ─┘   (mit Entkopplung)     └──► v_cq* ─┘
                                 ▲
 v_g,abc ─► Clarke ─► Park ──────┤        θ ◄── PLL (regelt v_q → 0)
 i_abc   ─► Clarke ─► Park ──────┘
```

### SRF-PLL

1. $v_{g,abc}$ → Clarke → Park mit geschätztem $\hat\theta$.
2. PI-Regler treibt $v_q \to 0$: $\hat\omega = \omega_0 + \mathrm{PI}(v_q)$.
3. $\hat\theta = \int \hat\omega\,dt$.

Wenn $v_q = 0$, zeigt die d-Achse exakt auf den Spannungszeiger.

---

## 6. Besonderheiten beim MMC (HVDC)

### 6.1 Arm-Größen zerlegen

Pro Phase $k \in \{a,b,c\}$: oberer Arm (Strom $i_{u,k}$, Spannung $v_{u,k}$),
unterer Arm ($i_{l,k}$, $v_{l,k}$), je Armdrossel $L_{arm}$.

$$
i_{s,k} = i_{u,k} - i_{l,k} \quad(\text{AC-Strom}),\qquad
i_{c,k} = \tfrac12(i_{u,k} + i_{l,k}) \quad(\text{Kreis-/Differenzstrom})
$$

$$
v_{s,k} = \tfrac12(v_{l,k} - v_{u,k}) \quad(\text{treibt AC}),\qquad
v_{c,k} = \tfrac{V_{dc}}{2} - \tfrac12(v_{u,k} + v_{l,k}) \quad(\text{treibt } i_c)
$$

Daraus entstehen **zwei entkoppelte Kreise**:

- **AC-Kreis:** $\left(L_T + \tfrac{L_{arm}}{2}\right)\frac{di_s}{dt} = v_s - v_g - \dots$
  → genau das Modell aus Abschnitt 5 mit $L = L_T + L_{arm}/2$.
  Der normale dq-Stromregler funktioniert 1:1.
- **Kreisstrom-Kreis:** $L_{arm}\frac{di_c}{dt} + R_{arm} i_c = v_c$.

### 6.2 Kreisstrom: dq mit −2ω (CCSC)

Der Kreisstrom $i_{c,k}$ enthält:
- **DC-Anteil** $I_{dc}/3$ (gleich in allen Phasen → **Nullsystem**!),
- eine **2. Harmonische im Gegensystem** (Folge a-c-b) durch die
  pulsierende Kondensatorenergie in den Armen.

Die 2. Harmonische erhöht Verluste und Kondensatorrippel. Der
Circulating Current Suppression Controller (Tu et al., 2011) transformiert
$i_{c,abc}$ mit Clarke und dann Park mit **Winkel $-2\theta$**:

$$
\underline i_{c,dq,-2} = \underline i_{c,\alpha\beta}\; e^{+j2\theta}
$$

In diesem Rahmen wird die −2ω-Komponente zu einer **Gleichgröße**;
zwei PI-Regler mit Sollwert 0 und Entkopplungsterm $\mp 2\omega L_{arm}$
unterdrücken sie. (Das Skript zeigt das numerisch.)

### 6.3 Nullsystem nutzen

- **DC-Anteil des Kreisstroms** = Nullsystem von $i_c$ → steuert den
  Energieaustausch DC-Seite ↔ Arme (Gesamtenergie-/Vdc-Regelung).
- **Nullsystemspannung** (z. B. 3. Harmonische, ~1/6 der Grundschwingung)
  im Modulator: erhöht den Aussteuerbereich um ≈15 %; sie verschwindet
  netzseitig, da der Trafo (Δ oder ungeerdeter Sternpunkt) kein
  Nullsystem überträgt.

### 6.4 Unsymmetrische Netze (Fehlerfall)

Bei unsymmetrischer Spannung gibt es Mit- **und** Gegensystem. Im
+ω-dq-Rahmen erscheint das Gegensystem als 2ω-Schwingung. Lösung:
**zwei** dq-Rahmen (+θ und −θ), Trennung per Notch-Filter bei 2ω oder
DDSRF/DSOGI-PLL, und getrennte Regler für Mit- und Gegensystemstrom
(Grid-Code-Anforderungen für Fault-Ride-Through).

---

## 7. Spickzettel

| Größe | Formel (amplitudeninvariant) |
|---|---|
| Clarke α | $\tfrac23(x_a - \tfrac12 x_b - \tfrac12 x_c)$ |
| Clarke β | $(x_b - x_c)/\sqrt3$ |
| Null | $(x_a+x_b+x_c)/3$ |
| Park | $x_d + jx_q = (x_\alpha + jx_\beta)e^{-j\theta}$ |
| Leistung | $P=\tfrac32(v_di_d+v_qi_q)$, $Q=\tfrac32(v_qi_d-v_di_q)$ |
| Kopplung | $+\omega L i_q$ (d), $-\omega L i_d$ (q) |
| MMC AC-Induktivität | $L_T + L_{arm}/2$ |
| CCSC-Rahmen | Winkel $-2\theta$ |

---

## 8. Übungen

1. Rechne von Hand: $v_a = 100\cos(\omega t + 30°)$ symmetrisch, $\theta=\omega t$.
   Wie groß sind $v_d, v_q$? *(Lösung: 86,6 und 50)*
2. Zeige, dass bei $v_q = 0$ und $i_q < 0$ (in dieser Konvention) der
   Umrichter Blindleistung **ins Netz liefert** ($Q>0$).
3. Was passiert mit $x_d, x_q$, wenn die PLL-Frequenz 1 Hz daneben liegt?
   *(Sie drehen langsam mit 1 Hz – genau das korrigiert die PLL.)*
4. Füge im Skript ein Gegensystem mit 10 % hinzu und beobachte die
   100-Hz-Schwingung in $v_d, v_q$.
5. Leite die Kopplungsterme des CCSC her: Rahmen mit $-2\omega$ →
   $L_{arm}\,d\underline i/dt = \dots + j2\omega L_{arm}\underline i$.

## Literatur

- A. Yazdani, R. Iravani: *Voltage-Sourced Converters in Power Systems*, Wiley 2010 (Kap. 4, 8 – Standardwerk für dq-Regelung).
- R. Teodorescu, M. Liserre, P. Rodríguez: *Grid Converters for Photovoltaic and Wind Power Systems*, Wiley 2011 (PLL, unsymmetrische Netze).
- Q. Tu, Z. Xu, L. Xu: „Reduced Switching-Frequency Modulation and Circulating Current Suppression for MMC“, IEEE Trans. Power Delivery, 2011.
- K. Sharifabadi et al.: *Design, Control and Application of Modular Multilevel Converters for HVDC Transmission Systems*, Wiley 2016.

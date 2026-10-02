"""Interactive app: planning and teaching tool for the resilient modulus (Mr) of SCBA-modified granular subbase.

Run:  streamlit run app/streamlit_app.py
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

import scba_core as sc

DOI = "https://doi.org/10.5281/zenodo.23101272"
REPO = "https://github.com/oscar1485/scba-resilient-modulus-rf"
LEVEL_LABEL = {0: "0 % (suelo/subbase sin ceniza)", 5: "5 % de SCBA", 10: "10 % de SCBA"}
COLORS = {0: "#48C9B0", 5: "#5DADE2", 10: "#AF7AC5"}

st.set_page_config(page_title="Mr–SCBA: simulador de módulo resiliente", page_icon="🧱", layout="wide")


@st.cache_resource(show_spinner="Entrenando modelos con los datos del estudio…")
def get_bundle():
    return sc.build_bundle()


B = get_bundle()


def scba_notice():
    st.info(
        "**Alcance.** El estudio ensayó **una probeta por contenido de ceniza (0, 5 y 10 %)** y 15 estados de esfuerzo. "
        "Por eso solo se pueden simular esos tres contenidos: entre ellos el modelo no está respaldado por datos. "
        "Los resultados son orientativos; no sustituyen un ensayo.")


# =====================================================================  PAGE 1
def page_plan():
    st.header("🧪 Plan de ensayo: valores sugeridos para el laboratorio")
    st.caption("Programa de 15 secuencias del estudio (coincide con AASHTO T 307 para base/subbase). "
               "Las cargas se calculan automáticamente a partir de los esfuerzos y del diámetro de la probeta.")
    scba_notice()
    c1, c2, c3 = st.columns(3)
    levels = c1.multiselect("Contenido de SCBA (%)", list(sc.SCBA_LEVELS), default=list(sc.SCBA_LEVELS))
    diam = c2.number_input("Diámetro de la probeta (mm)", 50.0, 200.0, 100.0, 5.0,
                           help="El estudio usó 100 mm (deducido de carga/esfuerzo en los datos).")
    gauge = c3.number_input("Longitud de medida de los LVDT (mm)", 20.0, 400.0, 200.0, 10.0,
                            help="Para estimar el desplazamiento recuperable que debe resolver el sensor.")
    if not levels:
        st.warning("Selecciona al menos un contenido de SCBA.")
        return
    tables = [sc.plan_table(B, int(cb), diam, gauge) for cb in levels]
    df = pd.concat(tables, ignore_index=True)
    st.dataframe(df, hide_index=True)
    st.caption(
        "**Mr medido** = mediana de las lecturas de esa secuencia en este estudio (una probeta por contenido; el rango P25–P75 "
        "refleja la variación entre lecturas de la misma probeta, no entre probetas). **k–θ** = ley potencial ajustada "
        "(Mr = k₁·θ^k₂). **εr** y **ΔL** se estiman con el Mr medido (k–θ si no hay dato válido). "
        "El Mr de la secuencia 15 con 10 % no se muestra: el LVDT axial marcó cero y esas lecturas se descartaron.")
    st.download_button("⬇️ Descargar tabla (CSV)", df.to_csv(index=False).encode("utf-8-sig"),
                       "plan_ensayo_scba.csv", "text/csv")

    st.subheader("Mapa del protocolo")
    s = B.sched
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.scatter(s["s3"], s["sd"], s=60, color="#2F5496", zorder=3)
    for _, r in s.iterrows():
        ax.annotate(int(r["seq"]), (r["s3"], r["sd"]), textcoords="offset points", xytext=(6, 4), fontsize=9)
    ax.set_xlabel("Esfuerzo de confinamiento σ3 (kPa)")
    ax.set_ylabel("Esfuerzo desviador cíclico σd (kPa)")
    ax.grid(True, linestyle="--", alpha=0.4)
    st.pyplot(fig)
    plt.close(fig)
    st.caption("Cada punto es una secuencia (número). Dentro de la envolvente de estos puntos el simulador es más confiable.")


# =====================================================================  PAGE 2
def page_sim():
    st.header("🔬 Simulador de un estado de esfuerzos")
    scba_notice()
    left, right = st.columns([1, 2])
    with left:
        cb = st.selectbox("Contenido de SCBA", list(sc.SCBA_LEVELS), format_func=lambda v: LEVEL_LABEL[v])
        mode = st.radio("Estado de esfuerzos", ["Secuencia del protocolo", "Personalizado"])
        diam = st.number_input("Diámetro de la probeta (mm)", 50.0, 200.0, 100.0, 5.0)
        gauge = st.number_input("Longitud de medida LVDT (mm)", 20.0, 400.0, 200.0, 10.0)
        if mode == "Secuencia del protocolo":
            seq = st.selectbox("Secuencia", list(range(1, 16)), index=7)
            row = B.sched[B.sched["seq"] == seq].iloc[0]
            s3, sd = float(row["s3"]), float(row["sd"])
            st.write(f"σ3 = **{s3:.1f} kPa**, σd = **{sd:.1f} kPa**")
        else:
            s3 = st.slider("σ3: confinamiento (kPa)", 20.7, 137.9, 68.9, 0.1)
            sd = st.slider("σd: desviador cíclico (kPa)", 18.6, 248.2, 125.0, 0.1)

    contact = float(sc.contact_from_cyclic(sd))
    smax = sd + contact
    theta = 3 * s3 + sd
    mr_rf = float(sc.predict_rf(B, cb, s3, sd)[0])
    mr_k = float(sc.predict_kth(B, cb, s3, sd)[0])

    with right:
        st.subheader("Cargas a aplicar")
        a, b_, c_, d_ = st.columns(4)
        a.metric("Carga cíclica (kN)", f"{float(sc.load_kN(sd, diam)):.3f}")
        b_.metric("Carga de contacto (kN)", f"{float(sc.load_kN(contact, diam)):.3f}")
        c_.metric("Carga máxima (kN)", f"{float(sc.load_kN(smax, diam)):.3f}")
        d_.metric("θ = 3σ3 + σd (kPa)", f"{theta:.0f}")

        st.subheader("Módulo resiliente estimado")
        e, f, g, h = st.columns(4)
        e.metric("Mr, Random Forest (MPa)", f"{mr_rf:.0f}", help=f"Error típico ≈ ±{sc.RMSE_RF_V2:.0f} MPa en estados no vistos")
        f.metric("Mr, ley k–θ (MPa)", f"{mr_k:.0f}", help=f"Error típico ≈ ±{sc.RMSE_KTH_V2:.0f} MPa en estados no vistos")
        eps = float(sc.resilient_strain_microstrain(sd, mr_k))
        g.metric("εr esperada (µε)", f"{eps:.0f}")
        h.metric("ΔL recuperable (mm)", f"{eps * 1e-6 * gauge:.3f}")

        inside = sc.in_envelope(B, s3, sd)
        near = sc.nearest_sequence(B, s3, sd)
        if inside:
            st.success(f"Estado dentro de la envolvente ensayada (secuencia más cercana: {near}).")
        else:
            st.warning("Estado **fuera** de la envolvente de los 15 estados ensayados: es extrapolación y no está validada.")
        if mode == "Secuencia del protocolo":
            st.caption("En los 15 estados del protocolo el Random Forest reproduce lo medido porque se entrenó con esos datos; "
                       "su utilidad está en los estados intermedios (modo *Personalizado*).")
        else:
            st.caption(f"En validación dejando una secuencia fuera, el error típico fue ≈ {sc.RMSE_RF_V2:.0f} MPa (RF) y "
                       f"≈ {sc.RMSE_KTH_V2:.0f} MPa (k–θ). El esfuerzo de contacto se toma como 10 % de σmáx.")

    t1, t2, t3 = st.tabs(["Mr vs θ", "Pulso de carga y respuesta (esquema)", "Comparar contenidos de SCBA"])

    with t1:
        d = B.data[B.data["% CBCA"] == cb]
        fig, ax = plt.subplots(figsize=(7.5, 4.5))
        ax.scatter(d["theta"], d[sc.TARGET], s=12, alpha=0.5, color=COLORS[cb], label="Lecturas medidas (depuradas)")
        th = np.linspace(80, 680, 100)
        a_, k2 = B.kth[cb]
        ax.plot(th, np.exp(a_ + k2 * np.log(th)), color="#333333", label="Ley k–θ ajustada")
        ax.scatter([theta], [mr_rf], marker="*", s=180, color="#E74C3C", zorder=5, label="RF (estado simulado)")
        ax.scatter([theta], [mr_k], marker="D", s=50, color="#333333", zorder=5, label="k–θ (estado simulado)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("θ = 3σ3 + σd (kPa)")
        ax.set_ylabel("Mr (MPa)")
        ax.grid(True, which="both", linestyle="--", alpha=0.3)
        ax.legend(fontsize=8)
        st.pyplot(fig)
        plt.close(fig)
        st.caption("En materiales granulares Mr crece con el esfuerzo volumétrico θ (endurecimiento por esfuerzo).")

    with t2:
        t = np.linspace(0, 2, 801)
        sig = sc.haversine_pulse(sd, t)
        eps_t = sig / (mr_k * 1000.0) * 1e6
        fig, ax1 = plt.subplots(figsize=(8, 4))
        ax1.plot(t, sig, color="#2F5496")
        ax1.set_xlabel("Tiempo (s)")
        ax1.set_ylabel("σd cíclico (kPa)", color="#2F5496")
        ax2 = ax1.twinx()
        ax2.plot(t, eps_t, color="#E74C3C", linestyle="--")
        ax2.set_ylabel("εr idealizada (µε)", color="#E74C3C")
        ax1.grid(True, linestyle="--", alpha=0.3)
        st.pyplot(fig)
        plt.close(fig)
        st.caption("Esquema didáctico: pulso haversine de 0,1 s seguido de 0,9 s de reposo (valores habituales en AASHTO T 307) y respuesta "
                   "elástica idealizada ε = σ/Mr. En un ensayo real aparece además deformación permanente acumulada y un pequeño retardo; "
                   "este modelo no las predice.")

    with t3:
        vals = [(c, float(sc.predict_rf(B, c, s3, sd)[0]), float(sc.predict_kth(B, c, s3, sd)[0])) for c in sc.SCBA_LEVELS]
        x = np.arange(3)
        fig, ax = plt.subplots(figsize=(6.5, 4))
        ax.bar(x - 0.2, [v[1] for v in vals], 0.4, label="Random Forest", color="#5DADE2", edgecolor="#333333")
        ax.bar(x + 0.2, [v[2] for v in vals], 0.4, label="k–θ", color="#F5D140", edgecolor="#333333")
        ax.set_xticks(x)
        ax.set_xticklabels([f"{c} %" for c in sc.SCBA_LEVELS])
        ax.set_ylabel("Mr (MPa)")
        ax.set_xlabel("Contenido de SCBA")
        ax.legend()
        ax.grid(True, axis="y", linestyle="--", alpha=0.3)
        st.pyplot(fig)
        plt.close(fig)
        st.warning("Con **una sola probeta por contenido**, las diferencias entre barras mezclan el efecto de la ceniza con la variabilidad "
                   "entre probetas. No se debe concluir un 'contenido óptimo' a partir de estas tres barras.")


# =====================================================================  PAGE 3
def page_learn():
    st.header("📚 Aprende: módulo resiliente axial")
    st.markdown("""
**¿Qué es el módulo resiliente (Mr)?** Es la rigidez de un material ante cargas repetidas pequeñas, calculada como el cociente
entre el esfuerzo desviador cíclico y la deformación axial **recuperable** (resiliente):

$$M_r = \\frac{\\Delta\\sigma_d}{\\varepsilon_r}$$

**¿Cómo es el ensayo?** Una probeta cilíndrica se confina con presión lateral σ3 y se le aplican pulsos axiales (σd). En cada pulso
se mide la deformación recuperable εr (Mr) y la deformación que no se recupera (deformación permanente). El protocolo del estudio
usa 15 secuencias que combinan 5 niveles de confinamiento con 3 niveles de σd.

**¿Por qué cambia Mr?** En materiales granulares Mr crece con el esfuerzo volumétrico θ = σ1 + 2σ3 = 3σ3 + σd
(modelo k–θ: Mr = k₁·θ^k₂). La ceniza de bagazo (SCBA) modifica la matriz y puede cambiar esa respuesta.
""")
    st.subheader("Calculadora: de esfuerzo y deformación a Mr")
    c1, c2, c3 = st.columns(3)
    sd = c1.number_input("σd cíclico (kPa)", 1.0, 500.0, 100.0, 1.0)
    eps = c2.number_input("εr recuperable (µε)", 10.0, 5000.0, 500.0, 10.0)
    c3.metric("Mr (MPa)", f"{sd / (eps * 1e-6) / 1000.0:.0f}")
    st.caption("1 µε = 10⁻⁶. Con 100 kPa y 500 µε resulta Mr = 200 MPa.")

    st.subheader("Parámetros de la ley k–θ ajustada (este estudio)")
    rows = [{"SCBA (%)": c, "k₁": float(np.exp(a)), "k₂": float(k2)} for c, (a, k2) in B.kth.items()]
    st.dataframe(pd.DataFrame(rows).round(3), hide_index=True)
    st.caption("Mr = k₁·θ^k₂ (Mr en MPa, θ en kPa). Un k₂ mayor significa más endurecimiento por esfuerzo; "
               "el valor del 10 % es muy distinto al de 0 % y 5 %.")

    st.subheader("Autoevaluación")
    q1 = st.radio("1. ¿Qué mide la deformación resiliente?", ["La deformación que se recupera al descargar",
                                                              "La deformación total acumulada", "El asentamiento por consolidación"],
                  index=None)
    if q1:
        st.write("✅ Correcto." if q1.startswith("La deformación que se recupera") else
                 "❌ Es la parte recuperable de cada pulso; la acumulada es la deformación permanente.")
    q2 = st.radio("2. Si aumenta el confinamiento σ3, en un material granular Mr normalmente…",
                  ["Aumenta", "Disminuye", "No cambia"], index=None)
    if q2:
        st.write("✅ Correcto: aumenta θ y con él la rigidez." if q2 == "Aumenta" else "❌ Revisa la ley k–θ: Mr crece con θ.")
    q3 = st.radio("3. ¿Por qué un modelo puede tener R² ≈ 0.99 y aun así generalizar mal?",
                  ["Porque las lecturas vecinas de una misma probeta se parecen y quedan en entrenamiento y prueba",
                   "Porque el R² siempre es engañoso", "Porque el modelo es demasiado simple"], index=None)
    if q3:
        st.write("✅ Correcto: es el problema de las lecturas repetidas (ver 'Acerca del modelo')." if q3.startswith("Porque las lecturas")
                 else "❌ Piensa en que 900 lecturas provienen de solo 3 probetas.")


# =====================================================================  PAGE 4
def page_about():
    st.header("ℹ️ Acerca del modelo")
    st.markdown(f"""
**Datos.** 900 lecturas = 3 probetas (0, 5 y 10 % de SCBA) × 15 secuencias × 20 lecturas. Se descartaron **{B.n_removed}** lecturas
inválidas del 10 % (sensor axial en cero; carga fuera de objetivo), quedando **{len(B.data)}**.

**Modelos.** (1) Ley potencial k–θ por contenido de SCBA. (2) Random Forest con entradas conocidas *antes* del ensayo
(% SCBA, σ3, σd, esfuerzo de contacto, θ); hiperparámetros: {sc.RF_PARAMS}. Ambos se entrenan al abrir la app con los datos del repositorio.

**Validación (artículo).** R² con distintos esquemas (875 lecturas):
""")
    st.dataframe(pd.DataFrame({
        "Esquema": ["V0 partición aleatoria (cota superior, no válida)", "V1 bloques fuera (media ± sd)",
                    "V2 secuencia fuera (principal)", "V3 contenido de SCBA fuera"],
        "Random Forest (M4)": ["0.98", "0.70 ± 0.30", "0.78", "0.49"],
        "Ley k–θ": ["0.87", "0.81 ± 0.13", "0.84", "0.47"]}), hide_index=True)
    st.markdown(f"""
**Limitaciones.**
- Una sola probeta por contenido de SCBA: no hay réplicas; no se puede afirmar generalización a otras probetas o materiales.
- Extrapolar a un contenido de SCBA no ensayado (V3) da R² ≈ 0.4–0.5: **por eso la app solo ofrece 0, 5 y 10 %**.
- El Random Forest no supera a la ley k–θ clásica; ambos tienen error típico de 30–40 MPa en estados no vistos.
- No predice deformación permanente ni vida a fatiga.

**Cómo citar.** Hernández-Atencia, Y., Gutiérrez-Portela, F., Diaz-Triana, O. A., Pulecio-Díaz, J., Forero-Muñoz, F. (2026).
*Code and data for: Resilient Modulus of SCBA-Modified Granular Subbases: Experimental and Random Forest Analysis* (v1.0.0). Zenodo. {DOI}

Código y datos: {REPO} · Licencias: código MIT, datos CC BY 4.0.
""")


# =====================================================================  NAV
PAGES = {"🧪 Plan de ensayo": page_plan, "🔬 Simulador": page_sim,
         "📚 Aprende": page_learn, "ℹ️ Acerca del modelo": page_about}
st.sidebar.title("Mr–SCBA")
choice = st.sidebar.radio("Navegación", list(PAGES))
st.sidebar.caption(f"Herramienta de apoyo para laboratorio y docencia.  \n[DOI]({DOI}) · [GitHub]({REPO})")
PAGES[choice]()

import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

st.set_page_config(
    page_title="PVT Team 16", 
    page_icon="🛢️",  
    layout="wide", 
    initial_sidebar_state="expanded", 
    menu_items={
        "About": "Dibuat oleh Kelompok 16 dengan anggota Said Muhamad Khazafa (12223009), Abraham Gwen Bramanti (12223027), Azzikra Selky Saefana Putra (12223031), Raisha Dewi Endriyanti (12223044), Daniel Syahputra Barus (12223074)."
    }
)

# Page selector
page = st.sidebar.selectbox("Select Page", ["Home", "Gas PVT Calculator", "Oil PVT Calculator", "Water PVT Calculator", "About"])

# Sidebar
st.sidebar.title("Input Data")
st.sidebar.header("General Data")
p = st.sidebar.number_input("Initial Reservoir Pressure (psia)", min_value=0.0, value=3000.0, step=100.0)
psc = st.sidebar.number_input("Standard Pressure (psia)", min_value=0.0, value=14.7, step=0.1)
Tf = st.sidebar.number_input("Reservoir Temperature (°F)", min_value=0.0, value=200.0, step=10.0)
Tr = Tf + 460
sg = st.sidebar.number_input("Gas Specific Gravity", min_value=0.0, value=0.7, step=0.1)

if page == "Home":
    st.title("🧪 Program PVT Calculator Sederhana")
    st.header("Created by Kelompok 16")
    st.markdown("---")

    st.success("Selamat datang di **Program PVT Calculator Sederhana**!")
    st.markdown(
        """
        Program ini dibuat sebagai bagian dari tugas Mini Project matakuliah TM2202 - Fluida Reservoir Praktikum ITB.

        Program ini dibangun dengan bahasa pemrograman Python dan framework Streamlit.

        Program ini diharapkan dapat membantu anda menghitung dan memvisualisasikan properti fluida reservoir dengan berbagai korelasi yang umum dalam Teknik Perminyakan.

        Semoga membantu!
    """
    )
    st.info("Untuk memulai, silakan pilih halaman pada navigasi sidebar disebelah pojok kiri atas")

elif page == "Gas PVT Calculator":
    st.header("**PVT Calculator - Gas Properties**")
    st.sidebar.header("Gas Data")
    co2 = st.sidebar.number_input("CO2 Impurities (%mole)", min_value=0.0, max_value=100.0, value=0.0, step=1.0)
    h2s = st.sidebar.number_input("H2S Impurities (%mole)", min_value=0.0, max_value=100.0, value=0.0, step=1.0)
    n2 = st.sidebar.number_input("N2 Impurities (%mole)", min_value=0.0, max_value=100.0, value=0.0, step=1.0)

    if co2 + h2s + n2 > 100:
        st.sidebar.error("Total mole percent of impurities must not exceed 100%.")

    st.sidebar.header("Select Correlation")
    gas_type = st.sidebar.selectbox("Gas Type", ["Natural Gas", "Condensate Gas"])
    crit_corr = st.sidebar.selectbox("Critical Properties Correlation", ["Sutton", "Standing"])
    z_corr = st.sidebar.selectbox("Z Factor Correlation", ["Dranchuk-Abou-Kassem", "Beggs & Brill", "Hall Yarborough"])

    def t_sutton(sg, co2, h2s, n2):
        t = 169.2 + 349.5 * sg - 74 * (sg ** 2)
        return t - 80 * co2 + 130 * h2s - 250 * n2

    def p_sutton(sg, co2, h2s, n2):
        p = 756.8 - 131 * sg - 3.6 * (sg ** 2)
        return p + 440 * co2 + 600 * h2s - 170 * n2

    def t_standing(sg, co2, h2s, n2, gas_type):
        if gas_type == "Natural Gas":
            t = 168 + 325 * sg - 12.5 * (sg ** 2)
        else:
            t = 187 + 330 * sg - 71.5 * (sg ** 2)
        return t - 80 * co2 + 130 * h2s - 250 * n2

    def p_standing(sg, co2, h2s, n2, gas_type):
        if gas_type == "Natural Gas":
            p = 677 + 15 * sg - 37.5 * (sg ** 2)
        else:
            p = 706 + 51.7 * sg - 11.1 * (sg ** 2)
        return p + 440 * co2 + 600 * h2s - 170 * n2

    def calculate_critical_props(sg, co2, h2s, n2, crit_corr, gas_type):
        if crit_corr == "Sutton":
            Tpc = t_sutton(sg, co2, h2s, n2)
            Ppc = p_sutton(sg, co2, h2s, n2)
        elif crit_corr == "Standing":
            Tpc = t_standing(sg, co2, h2s, n2, gas_type)
            Ppc = p_standing(sg, co2, h2s, n2, gas_type)
        return Tpc, Ppc

    def calculate_pr_tr(p, Tf, Tpc, Ppc):
        Tpr = (Tf + 460) / Tpc
        Ppr = p / Ppc
        return Tpr, Ppr

    def calculate_cg(tpr, z, ppr, ppc):
        from math import exp
        A1 = 0.3265
        A2 = -1.07
        A3 = -0.5339
        A4 = 0.01569
        A5 = -0.05165
        A6 = 0.5475
        A7 = -0.7361
        A8 = 0.1844
        A9 = 0.1056
        A10 = 0.6134
        A11 = 0.721
        ropr = 0.27 * ppr / (z * tpr)
        zperp = A1 + A2 / tpr + A3 / tpr ** 3 + A4 / tpr ** 4 + A5 / tpr ** 5
        zperp += 2 * ropr * (A6 + A7 / tpr + A8 / tpr ** 2)
        zperp -= 5 * ropr ** 4 * A9 * (A7 / tpr + A8 / tpr ** 2)
        zperp += 2 * A10 * ropr / tpr ** 3 * (1 + A11 * ropr ** 2 - A11 ** 2 * ropr ** 4) * exp(-A11 * ropr ** 2)
        Cpr = 1 / ppr - 0.27 / (z ** 2 * tpr) * (zperp / (1 + (ropr / z) * zperp))
        return Cpr / ppc * 1e5

    def gas_density(z, p, Tf, sg):
        R = 10.732
        Mair = 28.9647
        Tr = Tf + 460
        return (sg * p * Mair) / (z * R * Tr)

    def gas_viscosity(Tf, Dg, sg):
        from math import exp
        Mair = 28.9647 * sg
        K = (9.4 + 0.02 * Mair) * (Tf ** 1.5) / (209 + 19 * Mair + Tf)
        x = 3.5 + 986 / Tf + 0.01 * Mair
        Y = 2.4 - 0.2 * x
        return 1e-4 * K * exp(x * ((Dg / 62.4) ** Y))

    def z_dak(ppr, tpr):
        from math import exp
        def fRhoR(z):
            A11 = 0.721
            cons1 = 0.3265
            cons2 = -1.07
            cons3 = -0.5339
            cons4 = 0.01569
            cons5 = -0.05165
            cons6 = 0.5475
            cons7 = -0.7361
            cons8 = 0.1844
            cons9 = 0.1056
            cons10 = 0.6134

            Value1 = cons1 + cons2 / tpr + cons3 / tpr ** 3 + cons4 / tpr ** 4 + cons5 / tpr ** 5
            Value2 = 0.27 * ppr / tpr
            Value3 = cons6 + cons7 / tpr + cons8 / tpr ** 2
            Value4 = cons9 * (cons7 / tpr + cons8 / tpr ** 2)
            Value5 = cons10 / tpr ** 3
            rhoR = 0.27 * ppr / (z * tpr)
            return Value1 * rhoR - Value2 / rhoR + Value3 * rhoR ** 2 - Value4 * rhoR ** 5 + Value5 * (1 + A11 * rhoR ** 2) * rhoR ** 2 * exp(-A11 * rhoR ** 2) + 1

        def dfRho(z):
            return (fRhoR(z + 0.001) - fRhoR(z)) / 0.001

        z = 0.5
        for _ in range(150):
            fval = fRhoR(z)
            if abs(fval) < 1e-12:
                break
            z -= fval / dfRho(z)
        return z

    def z_beggs_brill(ppr, tpr):
        from math import log10, exp
        cons1 = 1.39 * ((tpr - 0.92) ** 0.5) - (0.36 * tpr) - 0.101
        cons2 = ((0.62 - 0.23 * tpr) * ppr) + (((0.066 / (tpr - 0.86)) - 0.037) * (ppr ** 2)) + ((0.32 / ((10 ** 9) ** (tpr - 1))) * (ppr ** 6))
        cons3 = (0.132 - (0.32 * log10(tpr)))
        cons4 = 10 ** (0.3016 - (0.49 * tpr) + (0.1825 * (tpr ** 2)))
        return cons1 + ((1 - cons1) / exp(cons2)) + (cons3 * (ppr ** cons4))

    def z_hall_yarborough(ppr, tpr):
        from math import exp
        t = 1 / tpr
        def FY(z):
            x1 = 0.06125 * ppr * t * exp(-1.2 * (1 - t) ** 2) / z
            x2 = -0.06125 * ppr * t * exp(-1.2 * ((1 - t) ** 2))
            x3 = 14.76 * t - 9.76 * (t ** 2) + 4.58 * (t ** 3)
            x4 = 90.7 * t - 242.2 * (t ** 2) + 42.4 * (t ** 3)
            x5 = 2.18 + 2.82 * t
            return x2 + ((x1 + (x1 ** 2) + (x1 ** 3) + (x1 ** 4)) / ((1 - x1) ** 3)) - (x3 * (x1 ** 2)) + (x4 * (x1 ** x5))

        def dFY(z):
            return (FY(z + 0.001) - FY(z)) / 0.001

        z = 0.5
        for _ in range(150):
            fval = FY(z)
            if abs(fval) < 1e-12:
                break
            z -= fval / dFY(z)
        return z

    def calculate_bg(z, Tf, p, psc=14.7, tsc=60):
        return 1000 * psc * z * (Tf + 460) / (5.61458 * (tsc + 460) * p)

    Tpc, Ppc = calculate_critical_props(sg, co2, h2s, n2, crit_corr, gas_type)
    Tpr, Ppr = calculate_pr_tr(p, Tf, Tpc, Ppc)

    if z_corr == "Dranchuk-Abou-Kassem":
        z = z_dak(Ppr, Tpr)
    elif z_corr == "Beggs & Brill":
        z = z_beggs_brill(Ppr, Tpr)
    elif z_corr == "Hall Yarborough":
        z = z_hall_yarborough(Ppr, Tpr)

    Bg = calculate_bg(z, Tf, p, psc)
    Cg = calculate_cg(Tpr, z, Ppr, Ppc)
    Dg = gas_density(z, p, Tf, sg)
    Vg = gas_viscosity(Tf, Dg, sg)

    st.write(f"**Tpc (°R)**: {Tpc:.2f}")
    st.write(f"**Ppc (psia)**: {Ppc:.2f}")
    st.write(f"**Tpr**: {Tpr:.4f}")
    st.write(f"**Ppr**: {Ppr:.4f}")
    st.success(f"Z-Factor: {z:.5f}")
    st.success(f"Bg: {Bg:.5f} (RB/scf)")
    st.success(f"Cg x 1e-5: {Cg:.5f} (1/MMpsi)")
    st.success(f"Gas Density: {Dg:.2f} lbm/cf")
    st.success(f"Gas Viscosity: {Vg:.5f} cP")

    
    summary = f"""Gas Properties Summary

    Tpc (°R): {Tpc:.2f}
    Ppc (psia): {Ppc:.2f}
    Tpr: {Tpr:.4f}
    Ppr: {Ppr:.4f}

    Z-Factor: {z:.5f}
    Bg: {Bg:.5f} RB/scf
    Cg x 1e-5: {Cg:.5f} (1/MMpsi)
    Gas Density: {Dg:.2f} lbm/ft³
    Gas Viscosity: {Vg:.5f} cP
    """

    summary_bytes = summary.encode("utf-8")


    st.download_button(
        label="Download Summary TXT",
        data=summary_bytes,
        file_name="gas_properties_summary.txt",
        mime="text/plain"
    )

    st.subheader("Gas Properties Table")
    pressures_table = np.linspace(p, psc, 37)
    z_values, bg_values, cg_values, dg_values, vg_values = [], [], [], [], []

    for pi in pressures_table:
        Tpr_i, Ppr_i = calculate_pr_tr(pi, Tf, Tpc, Ppc)
        if z_corr == "Dranchuk-Abou-Kassem":
            z_i = z_dak(Ppr_i, Tpr_i)
        elif z_corr == "Beggs & Brill":
            z_i = z_beggs_brill(Ppr_i, Tpr_i)
        elif z_corr == "Hall Yarborough":
            z_i = z_hall_yarborough(Ppr_i, Tpr_i)
        z_values.append(z_i)
        bg_values.append(calculate_bg(z_i, Tf, pi, psc))
        cg_values.append(calculate_cg(Tpr_i, z_i, Ppr_i, Ppc))
        dg_values.append(gas_density(z_i, pi, Tf, sg))
        vg_values.append(gas_viscosity(Tf, dg_values[-1], sg))

    df_gas = pd.DataFrame({
        "Pressure (psia)": pressures_table,
        "Z-Factor": z_values,
        "Bg (RB/scf)": bg_values,
        "Cg x 1e-5 (1/MMpsi)": cg_values,
        "Gas Density (lbm/ft³)": dg_values,
        "Gas Viscosity (cP)": vg_values
    })

    st.dataframe(df_gas)
    csv = df_gas.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="Download Table CSV",
        data=csv,
        file_name="gas_properties_table.csv",
        mime="text/csv",
    )

    st.subheader("Gas Properties Chart")
    chart_option = st.selectbox("Select Chart", ["Z-Factor", "Bg", "Cg", "Gas Density", "Gas Viscosity"])

    pressures = np.linspace(psc, p, 100)
    fig, ax = plt.subplots()

    if chart_option == "Z-Factor":
        z_list = []
        for pi in pressures:
            Tpr_i, Ppr_i = calculate_pr_tr(pi, Tf, Tpc, Ppc)
            if z_corr == "Dranchuk-Abou-Kassem":
                z_i = z_dak(Ppr_i, Tpr_i)
            elif z_corr == "Beggs & Brill":
                z_i = z_beggs_brill(Ppr_i, Tpr_i)
            elif z_corr == "Hall Yarborough":
                z_i = z_hall_yarborough(Ppr_i, Tpr_i)
            z_list.append(z_i)
        ax.plot(pressures, z_list)
        ax.set_ylabel("Z-Factor")
        ax.set_title("Z-Factor vs Pressure")

    elif chart_option == "Bg":
        bg_list = []
        for pi in pressures:
            Tpr_i, Ppr_i = calculate_pr_tr(pi, Tf, Tpc, Ppc)
            if z_corr == "Dranchuk-Abou-Kassem":
                z_i = z_dak(Ppr_i, Tpr_i)
            elif z_corr == "Beggs & Brill":
                z_i = z_beggs_brill(Ppr_i, Tpr_i)
            elif z_corr == "Hall Yarborough":
                z_i = z_hall_yarborough(Ppr_i, Tpr_i)
            bg_list.append(calculate_bg(z_i, Tf, pi, psc))
        ax.plot(pressures, bg_list)
        ax.set_ylabel("Bg (RB/scf)")
        ax.set_title("Gas Formation Volume Factor (Bg) vs Pressure")

    elif chart_option == "Cg":
        cg_list = []
        for pi in pressures:
            Tpr_i, Ppr_i = calculate_pr_tr(pi, Tf, Tpc, Ppc)
            if z_corr == "Dranchuk-Abou-Kassem":
                z_i = z_dak(Ppr_i, Tpr_i)
            elif z_corr == "Beggs & Brill":
                z_i = z_beggs_brill(Ppr_i, Tpr_i)
            elif z_corr == "Hall Yarborough":
                z_i = z_hall_yarborough(Ppr_i, Tpr_i)
            cg_list.append(calculate_cg(Tpr_i, z_i, Ppr_i, Ppc))
        ax.plot(pressures, cg_list)
        ax.set_ylabel("Cg x 1e-5 (1/MMpsi)")
        ax.set_title("Gas Compressibility (Cg) vs Pressure")

    elif chart_option == "Gas Density":
        dg_list = []
        for pi in pressures:
            Tpr_i, Ppr_i = calculate_pr_tr(pi, Tf, Tpc, Ppc)
            if z_corr == "Dranchuk-Abou-Kassem":
                z_i = z_dak(Ppr_i, Tpr_i)
            elif z_corr == "Beggs & Brill":
                z_i = z_beggs_brill(Ppr_i, Tpr_i)
            elif z_corr == "Hall Yarborough":
                z_i = z_hall_yarborough(Ppr_i, Tpr_i)
            dg_list.append(gas_density(z_i, pi, Tf, sg))
        ax.plot(pressures, dg_list)
        ax.set_ylabel("Gas Density (lbm/ft³)")
        ax.set_title("Gas Density vs Pressure")

    elif chart_option == "Gas Viscosity":
        mu_list = []
        for pi in pressures:
            Tpr_i, Ppr_i = calculate_pr_tr(pi, Tf, Tpc, Ppc)
            if z_corr == "Dranchuk-Abou-Kassem":
                z_i = z_dak(Ppr_i, Tpr_i)
            elif z_corr == "Beggs & Brill":
                z_i = z_beggs_brill(Ppr_i, Tpr_i)
            elif z_corr == "Hall Yarborough":
                z_i = z_hall_yarborough(Ppr_i, Tpr_i)
            dg_i = gas_density(z_i, pi, Tf, sg)
            mu_list.append(gas_viscosity(Tf, dg_i, sg))
        ax.plot(pressures, mu_list)
        ax.set_ylabel("Viscosity (cP)")
        ax.set_title("Gas Viscosity vs Pressure")

    ax.set_xlabel("Pressure (psia)")
    ax.set_xlim(0, None)
    st.pyplot(fig)

elif page == "Water PVT Calculator":
    st.header("**PVT Calculator - Oil Properties**")
    st.sidebar.header("Water Data")
    ws = st.sidebar.number_input("TDS (wt%)", min_value=0.0, max_value=100.0, value=5.0, step=1.0)
    sat_cond = st.sidebar.selectbox("Saturation Condition", ["Gas-Free Water", "Gas-Saturated Water"])

    def bw(Tf, p, sat_cond):
        if sat_cond == 'Gas-Free Water':
            a1_A1, a2_A1, a3_A1 = 0.9947, 5.8e-6, 1.02e-6
            a1_A2, a2_A2, a3_A2 = -4.228e-6, 1.8376e-8, -6.77e-11
            a1_A3, a2_A3, a3_A3 = 1.3e-10, -1.3855e-12, 4.285e-15
        else:
            a1_A1, a2_A1, a3_A1 = 0.9911, 6.35e-5, 8.5e-7
            a1_A2, a2_A2, a3_A2 = -1.093e-6, -3.497e-9, 4.57e-12
            a1_A3, a2_A3, a3_A3 = -5.0e-11, 6.429e-13, -1.43e-15

        A1 = a1_A1 + a2_A1 * Tf + a3_A1 * Tf ** 2
        A2 = a1_A2 + a2_A2 * Tf + a3_A2 * Tf ** 2
        A3 = a1_A3 + a2_A3 * Tf + a3_A3 * Tf ** 2

        return A1 + A2 * p + A3 * p ** 2, A1, A2, A3
    
    def water_dens(ws, Bw):
        rho_w = 62.368 + 0.438603 * ws + 1.6007e-3 * ws ** 2
        return rho_w / Bw

    def water_visc(ws, Tf, p):
        a = 4.518e-2 + 9.313e-7 * ws - 3.93e-12 * ws ** 2
        b = 70.634 + 9.576e-10 * ws ** 2
        wvis = a + b / Tf
        return ((1 + 3.5e-2 * p ** 2 * (Tf - 40)) * wvis) * 1e-8

    def rsw(Tf, ws, p):
        A0, A1, A2, A3 = 8.15839, -0.0612265, 0.000191663, -2.1654e-7
        B0, B1, B2, B3 = 1.01021e-2, -7.44241e-5, 3.05553e-7, -2.94883e-10
        C0, C1, C2, C3, C4 = -9.02505, 0.130237, -8.53425e-4, 2.34122e-6, -2.37049e-9

        a = A0 + A1 * Tf + A2 * Tf**2 + A3 * Tf**3
        b = B0 + B1 * Tf + B2 * Tf**2 + B3 * Tf**3
        C = (C0 + C1 * Tf + C2 * Tf**2 + C3 * Tf**3 + C4 * Tf**4) * 1e-7

        Rsw = (a + b * p + C * p**2) * 1e-3
        x = -0.0840655 * ws * (Tf ** -0.285854)
        Y = 10 ** x
        return Rsw * Y

    def cw(p, Tf, Rsw, ws, sat_cond):
        a = 3.8546 - 0.000134 * p
        b = -0.01052 + 4.77e-7 * p
        C = 3.9267e-5 - 8.8e-10 * p
        Cwt = (a + b * Tf + C * Tf**2) / 1e6 * 1e5

        if sat_cond == "Gas-Free Water":
            return Cwt
        else:
            return Cwt * (1 + 8.9e-3 * Rsw)

    Bw, A1, A2, A3 = bw(Tf, p, sat_cond)
    st.success(f"Bw = {Bw:.5f} (RB/STB)")

    rho_w = water_dens(ws, Bw)
    st.success(f"Water Density = {rho_w:.2f} (lbm/cf)")

    mu_w = water_visc(ws, Tf, p)
    st.success(f"Water Visc. = {mu_w:.5f} (cp)")

    Rsw = rsw(Tf, ws, p)
    st.success(f"Rsw = {Rsw:.5f} (scf/STB)")

    Cw = cw(p, Tf, Rsw, ws, sat_cond)
    st.success(f"Cw = {Cw:.5f} (1/MMpsi)")

    water_summary_txt = f"""Water Properties Summary

    Bw (RB/STB): {Bw:.5f}
    Water Density (lbm/cf): {rho_w:.2f}
    Water Viscosity (cP): {mu_w:.5f}
    Rsw (scf/STB): {Rsw:.5f}
    Cw (1/MMpsi): {Cw:.5f}
    """

    water_summary_bytes = water_summary_txt.encode("utf-8")

    st.download_button(
        label="Download Summary TXT",
        data=water_summary_bytes,
        file_name="water_properties_summary.txt",
        mime="text/plain"
    )

    st.subheader("Water Properties Table")
    pressures_table = np.linspace(p, psc, 37)
    Bw_list, rho_w_list, mu_w_list, Rsw_list, Cw_list = [], [], [], [], []

    for pi in pressures_table:
        Bw_i, _, _, _ = bw(Tf, pi, sat_cond)
        rho_w_i = water_dens(ws, Bw_i)
        mu_w_i = water_visc(ws, Tf, pi)
        Rsw_i = rsw(Tf, ws, pi)
        Cw_i = cw(pi, Tf, Rsw_i, ws, sat_cond)

        Bw_list.append(Bw_i)
        rho_w_list.append(rho_w_i)
        mu_w_list.append(mu_w_i)
        Rsw_list.append(Rsw_i)
        Cw_list.append(Cw_i)

    df_water = pd.DataFrame({
        "Pressure (psia)": pressures_table,
        "Bw (RB/STB)": Bw_list,
        "Water Density (lbm/ft³)": rho_w_list,
        "Water Viscosity (cP)": mu_w_list,
        "Rsw (scf/STB)": Rsw_list,
        "Cw (1/MMpsi)": Cw_list
    })

    st.dataframe(df_water)
    csv = df_water.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="Download Table CSV",
        data=csv,
        file_name="water_properties_table.csv",
        mime="text/csv",
    )

    st.subheader("Water Properties Chart")
    chart_option = st.selectbox("Select Chart", ["Bw", "Water Dens", "Water Visc", "Rsw", "Cw"])

    pressures = np.linspace(psc, p, 100)
    fig, ax = plt.subplots()

    if chart_option == "Bw":
        _, A1, A2, A3 = bw(Tr, p, sat_cond)
        Bws = A1 + A2 * pressures + A3 * pressures**2
        ax.plot(pressures, Bws)
        ax.set_ylabel("Bw (RB/STB)")
        ax.set_title("Water FVF vs Pressure")
    elif chart_option == "Water Dens":
        Bws = [bw(Tf, pi, sat_cond)[0] for pi in pressures]
        densities = [water_dens(ws, bwi) for bwi in Bws]
        ax.plot(pressures, densities)
        ax.set_ylabel("Water Density (lb/ft³)")
        ax.set_title("Water Density vs Pressure")
    elif chart_option == "Water Visc":
        viscosities = [water_visc(ws, Tf, pi) for pi in pressures]
        ax.plot(pressures, viscosities)
        ax.set_ylabel("Water Viscosity (cp)")
        ax.set_title("Water Viscosity vs Pressure")
    elif chart_option == "Rsw":
        rsws = [rsw(Tf, ws, pi) for pi in pressures]
        ax.plot(pressures, rsws)
        ax.set_ylabel("Rsw (SCF/STB)")
        ax.set_title("Gas Solubility (Rsw) vs Pressure")
    elif chart_option == "Cw":
        rsws = [rsw(Tf, ws, pi) for pi in pressures]
        cws = [cw(pi, Tf, rsws[i], ws, sat_cond) for i, pi in enumerate(pressures)]
        ax.plot(pressures, cws)
        ax.set_ylabel("Cw (1/MMpsi)")
        ax.set_title("Water Compressibility (Cw) vs Pressure")

    ax.set_xlabel("Pressure (psia)")
    ax.set_xlim(0, None)
    st.pyplot(fig)

elif page == "Oil PVT Calculator":
    st.write("**PVT Calculator - Oil Properties**")
    st.sidebar.header("Oil Data")
    api = st.sidebar.number_input("Oil API Gravity (°API)", min_value=0.0, value=35.0, step=1.0)
    sgoil = st.sidebar.number_input("Oil Specific Gravity", min_value=0.0, value=0.85, step=0.1)
    pb = st.sidebar.number_input("Bubble Point Pressure (psia)", min_value=0.0, value=2500.0, step=100.0)
    psep = st.sidebar.number_input("Separator Pressure (psia)", min_value=0.0, value=100.0, step=10.0)
    tsep = st.sidebar.number_input("Separator Temperature (°F)", min_value=0.0, value=100.0, step=10.0)
    st.sidebar.header("Select Correlation")
    rs_method = st.sidebar.selectbox("Rs Correlation", ["Standing", "Vazquez-Beggs"])
    rho_method = st.sidebar.selectbox("Oil Density Correlation", ["Standard", "Standing"])

    condition = "Saturated" if p < pb else "Undersaturated"
    st.info(f"Reservoir condition: {condition}")

    def rs_standing(sg, p, t, api, pb):
        x = 0.0125 * api - 0.00091 * t
        target_p = p if p < pb else pb
        return sg * (((target_p / 18.2) + 1.4) * (10 ** x)) ** 1.2048

    def rs_vazquez_beggs(sg, p, t, api, tsep, psep, pb):
        sggs = sg * (1 + 5.912e-5 * api * tsep * np.log10(psep / 114.7))
        if api <= 30:
            C1, C2, C3 = 0.0362, 1.0937, 25.724
        else:
            C1, C2, C3 = 0.0178, 1.187, 23.931
        target_p = p if p < pb else pb
        return C1 * sggs * (target_p ** C2) * np.exp(C3 * (api / (t + 460)))

    def coil(rs, t, api, sg, p, pb, tsep, psep):
        sggs = sg * (1 + 5.912e-5 * api * tsep * np.log10(psep / 114.7))
        if p >= pb:
            return (-1433 + 5 * rs + 17.2 * t - 1180 * sggs + 12.61 * api) / (1e5 * p)
        else:
            a = -7.573 - 1.45 * np.log(p) - 0.383 * np.log(pb) + 1.402 * np.log(t + 460) + 0.256 * np.log(api) + 0.449 * np.log(rs)
            return np.exp(a)

    def boil(rs, sg, sgoil, t, pb, p, co):
        if p < pb:
            return 0.9759 + 0.00012 * ((rs * ((sg / sgoil) ** 0.5) + 1.25 * t) ** 1.2)
        else:
            bob = 0.9759 + 0.00012 * ((rs * ((sg / sgoil) ** 0.5) + 1.25 * t) ** 1.2)
            return bob * np.exp(-co * (p - pb))

    def oil_density_standar(sgoil, sg, Rs, Bo):
        return (62.4 * sgoil + 0.0136 * Rs * sg) / Bo

    def oil_density_standing(co, sgoil, sg, Rs, t, p, pb):
        if p < pb:
            return (62.4 * sgoil + 0.0136 * Rs * sg) / (0.972 + 0.000147 * ((Rs * ((sg / sgoil) ** 0.5) + 1.25 * t) ** 1.175))
        else:
            roob = (62.4 * sgoil + 0.0136 * Rs * sg) / (0.972 + 0.000147 * ((Rs * ((sg / sgoil) ** 0.5) + 1.25 * t) ** 1.175))
            return roob * np.exp(co * (p - pb))

    def dead_oil_viscosity(api, t):
        Zb = 3.0324 - 0.02023 * api
        Y = 10 ** Zb
        x = Y * (t ** -1.163)
        return (10 ** x) - 1

    def oil_viscosity_sat(Rs, api, t):
        a = 10.715 * ((Rs + 100) ** -0.515)
        b = 5.44 * ((Rs + 150) ** -0.338)
        return a * (dead_oil_viscosity(api, t) ** b)

    def oil_viscosity_unsat(p, pb, Rs, api, t):
        a = -3.9e-5 * p - 5
        m = 2.6 * (p ** 1.187) * (10 ** a)
        return oil_viscosity_sat(Rs, api, t) * ((p / pb) ** m)

    if rs_method == "Standing":
        Rs = rs_standing(sg, p, Tf, api, pb)
    elif rs_method == "Vazquez-Beggs":
        Rs = rs_vazquez_beggs(sg, p, Tf, api, tsep, psep, pb)

    st.success(f"Rs = {Rs/1000:.5f} (Mscf/STB)")

    co = coil(Rs, Tf, api, sg, p, pb, tsep, psep)
    st.success(f"Co x 1e-5 = {co * 1e5:.5f} (1/MMpsi)")

    Bo = boil(Rs, sg, sgoil, Tf, pb, p, co)
    st.success(f"Bo = {Bo:.5f} (RB/STB)")

    if rho_method == "Standard":
        rho_o = oil_density_standar(sgoil, sg, Rs, Bo)
    elif rho_method == "Standing":
        rho_o = oil_density_standing(co, sgoil, sg, Rs, Tf, p, pb)

    st.success(f"Oil Density = {rho_o:.2f} (lbm/cf)")

    if p < pb:
        mu_o = oil_viscosity_sat(Rs, api, Tf)
    else:
        mu_o = oil_viscosity_unsat(p, pb, Rs, api, Tf)
    st.success(f"Oil Viscosity = {mu_o:.3f} (cP)")

    oil_summary_txt = f"""Oil Properties Summary

    Reservoir Condition: {condition}

    Rs (Mscf/STB): {Rs / 1000:.5f}
    Co x 1e-5 (1/MMpsi): {co * 1e5:.5f}
    Bo (RB/STB): {Bo:.5f}
    Oil Density (lbm/cf): {rho_o:.2f}
    Oil Viscosity (cP): {mu_o:.3f}
    """

    st.download_button(
        label="Download Summary TXT",
        data=oil_summary_txt.encode("utf-8"),
        file_name="oil_properties_summary.txt",
        mime="text/plain"
    )

    st.subheader("Oil Properties Table")
    pressures_table = np.linspace(p, psc, 37)
    rs_values, co_values, bo_values, rhoo_values, muo_values = [], [], [], [], []

    for pi in pressures_table:
        if rs_method == "Standing":
            rs_i = rs_standing(sg, pi, Tf, api, pb)
        elif rs_method == "Vazquez-Beggs":
            rs_i = rs_vazquez_beggs(sg, pi, Tf, api, tsep, psep, pb)

        co_i = coil(rs_i, Tf, api, sg, pi, pb, tsep, psep)
        bo_i = boil(rs_i, sg, sgoil, Tf, pb, pi, co_i)

        if rho_method == "Standard":
            rho_i = oil_density_standar(sgoil, sg, rs_i, bo_i)
        elif rho_method == "Standing":
            rho_i = oil_density_standing(co_i, sgoil, sg, rs_i, Tf, pi, pb)

        if pi < pb:
            mu_i = oil_viscosity_sat(rs_i, api, Tf)
        else:
            mu_i = oil_viscosity_unsat(pi, pb, rs_i, api, Tf)

        rs_values.append(rs_i)
        co_values.append(co_i * 1e5)
        bo_values.append(bo_i)
        rhoo_values.append(rho_i)
        muo_values.append(mu_i)

    df_oil = pd.DataFrame({
        "Pressure (psia)": pressures_table,
        "Rs (scf/STB)": rs_values,
        "Co x 1e-5 (1/MMpsi)": co_values,
        "Bo (RB/STB)": bo_values,
        "Oil Density (lbm/ft³)": rhoo_values,
        "Oil Viscosity (cP)": muo_values
    })

    st.dataframe(df_oil)

    csv = df_oil.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="Download Table CSV",
        data=csv,
        file_name="oil_properties_table.csv",
        mime="text/csv",
    )



    st.subheader("Oil Properties Chart")
    chart_option = st.selectbox("Select Chart", ["Rs", "Bo", "Co", "Oil Density", "Oil Viscosity"])
    pressures = np.linspace(psc, p, 100)

    fig, ax = plt.subplots()
    if chart_option == "Rs":
        rs_curve = [rs_standing(sg, pi, Tf, api, pb) if rs_method == "Standing" else rs_vazquez_beggs(sg, pi, Tf, api, tsep, psep, pb) for pi in pressures]
        ax.plot(pressures, [r/1000 for r in rs_curve])
        ax.set_ylabel("Rs (Mscf/STB)")
        ax.set_title("Solution GOR (Rs) vs Pressure")
    elif chart_option == "Bo":
        rs_curve = [rs_standing(sg, pi, Tf, api, pb) if rs_method == "Standing" else rs_vazquez_beggs(sg, pi, Tf, api, tsep, psep, pb) for pi in pressures]
        co_curve = [coil(rs_curve[i], Tf, api, sg, pi, pb, tsep, psep) for i, pi in enumerate(pressures)]
        bo_curve = [boil(rs_curve[i], sg, sgoil, Tf, pb, pi, co_curve[i]) for i, pi in enumerate(pressures)]
        ax.plot(pressures, bo_curve)
        ax.set_ylabel("Bo (RB/STB)")
        ax.set_title("Oil FVF (Bo) vs Pressure")
    elif chart_option == "Co":
        rs_curve = [rs_standing(sg, pi, Tf, api, pb) if rs_method == "Standing" else rs_vazquez_beggs(sg, pi, Tf, api, tsep, psep, pb) for pi in pressures]
        co_curve = [coil(rs_curve[i], Tf, api, sg, pi, pb, tsep, psep) * 1e5 for i, pi in enumerate(pressures)]
        ax.plot(pressures, co_curve)
        ax.set_ylabel("Co x 1e-5 (1/MMpsi)")
        ax.set_title("Oil Compressibility (Co) vs Pressure")
    elif chart_option == "Oil Density":
        rs_curve = [rs_standing(sg, pi, Tf, api, pb) if rs_method == "Standing" else rs_vazquez_beggs(sg, pi, Tf, api, tsep, psep, pb) for pi in pressures]
        co_curve = [coil(rs_curve[i], Tf, api, sg, pi, pb, tsep, psep) for i, pi in enumerate(pressures)]
        if rho_method == "Standard":
            bo_curve = [boil(rs_curve[i], sg, sgoil, Tf, pb, pi, co_curve[i]) for i in range(len(pressures))]
            rho_curve = [oil_density_standar(sgoil, sg, rs_curve[i], bo_curve[i]) for i in range(len(pressures))]
        else:
            rho_curve = [oil_density_standing(co_curve[i], sgoil, sg, rs_curve[i], Tf, pi, pb) for i, pi in enumerate(pressures)]
        ax.plot(pressures, rho_curve)
        ax.set_ylabel("Oil Density (lbm/cf)")
        ax.set_title("Oil Density vs Pressure")
    elif chart_option == "Oil Viscosity":
        rs_curve = [rs_standing(sg, pi, Tf, api, pb) if rs_method == "Standing" else rs_vazquez_beggs(sg, pi, Tf, api, tsep, psep, pb) for pi in pressures]
        mu_curve = [oil_viscosity_sat(rs_curve[i], api, Tf) if pi < pb else oil_viscosity_unsat(pi, pb, rs_curve[i], api, Tf) for i, pi in enumerate(pressures)]
        ax.plot(pressures, mu_curve)
        ax.set_ylabel("Oil Viscosity (cP)")
        ax.set_title("Oil Viscosity vs Pressure")

    ax.set_xlabel("Pressure (psia)")
    ax.set_xlim(0, None)
    st.pyplot(fig)

elif page == "About":
    st.title("About Team 16")
    st.markdown(
        """
        Dosen pengampu: **Zuher Syihab, S.T., Ph.D.**

        Program ini dikembangkan oleh Kelompok 16 dengan anggota:

        1. Said Muhamad Khazafa (12223009)

        2. Abraham Gwen Bramanti (12223027)

        3. Azzikra Selky Saefana Putra (12223031)

        4. Raisha Dewi Endriyanti (12223044)

        5. Daniel Syahputra Barus (12223074)

        """
    )
    st.info("Terakhir diupdate 14 Juni 2025")
    st.caption("© 2025 Kelompok 16 - Teknik Perminyakan ITB")


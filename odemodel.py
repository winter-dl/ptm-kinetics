import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import matplotlib.cm as cm


def sim(variables, t, params):
    A = variables[0]
    B = variables[1]
    Ap = variables[2]
    AB = variables[3]
    ApB = variables[4]

    kon1 = params[0]
    koff1 = params[1]
    kon2 = params[2]
    koff2 = params[3]
    cat_kin1 = params[4]
    cat_pp1 = params[5]
    cat_kin2 = params[6]
    cat_pp2 = params[7]

    dAdt = -(kon1*A*B) + (koff1*AB) - (cat_kin1*A) + (cat_pp1*Ap)
    dBdt = -(kon1*A*B) + (koff1*AB) - (kon2*Ap*B) + (koff2*ApB)
    dApdt = -(kon2*Ap*B) + (koff2*ApB) + (cat_kin1*A) - (cat_pp1*Ap)
    dABdt = (kon1*A*B) - (koff1*AB) - (cat_kin2*AB) + (cat_pp2*ApB)
    dApBdt = (kon2*Ap*B) - (koff2*ApB) + (cat_kin2*AB) - (cat_pp2*ApB)

    return([dAdt, dBdt, dApdt, dABdt, dApBdt])


def ode_switch_model(y0: dict, kinetic_params: dict, eq_duration=200, input1_duration=600, input2_duration=600, inputs=2, pp_background="without_pp_background"):
    if inputs != 1 and inputs != 2:
        raise ValueError("inputs must be either 1 or 2.")
    if not pp_background in ["with_pp_background", "without_pp_background"]:
        raise ValueError(
            "pp_background must be either 'with_pp_background' or 'without_pp_background'.")
    if inputs == 2 and pp_background == "with_pp_background":
        raise ValueError(
            "It is not possible to set inputs to 2 and pp_background to 'with_pp_background' simultaneously.")

    eq_params = [kinetic_params["kon1"], kinetic_params["koff1"],
                 kinetic_params["kon2"], kinetic_params["koff2"], 0, 0, 0, 0]
    if pp_background == "with_pp_background":
        eq_params[4] = kinetic_params['cat_pp1']
        eq_params[6] = kinetic_params['cat_pp2']
    input1_params = eq_params[0:4]+[kinetic_params["cat_kin1"]] + \
        [eq_params[5]]+[kinetic_params["cat_kin2"]]+[eq_params[7]]
    input2_params = input1_params[0:5]+[kinetic_params["cat_pp1"]
                                        ]+[input1_params[6]]+[kinetic_params["cat_pp2"]]

    t_eq = np.linspace(0, eq_duration, num=eq_duration*60)
    t_input1 = np.linspace(eq_duration, eq_duration +
                           input1_duration, num=input1_duration*60)
    t_input2 = np.linspace(eq_duration + input1_duration, eq_duration +
                           input1_duration + input2_duration, num=input2_duration*60)

    # Calculate the equilibration phase

    y0_list = [y0["A"], y0["B"], y0["Ap"], y0["AB"], y0["ApB"]]
    y_eq = odeint(sim, y0_list, t_eq, args=(eq_params,))

    # Update variables

    y_after_eq = [y_eq[:, 0][-1], y_eq[:, 1][-1], y_eq[:, 2]
                  [-1], y_eq[:, 3][-1], y_eq[:, 4][-1]]

    # Calculate the input phase

    y_input1 = odeint(sim, y_after_eq, t_input1, args=(input1_params,))

    # Update variables again

    y_after_input1 = [y_input1[:, 0][-1], y_input1[:, 1][-1], y_input1[:, 2]
                      [-1], y_input1[:, 3][-1], y_input1[:, 4][-1]]

    # Calculate the second input phase

    y_input2 = odeint(sim, y_after_input1, t_input2, args=(input2_params,))

    t = np.concatenate((t_eq, t_input1, t_input2))
    A_values = np.concatenate((y_eq[:, 0], y_input1[:, 0], y_input2[:, 0]))
    B_values = np.concatenate((y_eq[:, 1], y_input1[:, 1], y_input2[:, 1]))
    Ap_values = np.concatenate((y_eq[:, 2], y_input1[:, 2], y_input2[:, 2]))
    AB_values = np.concatenate((y_eq[:, 3], y_input1[:, 3], y_input2[:, 3]))
    ApB_values = np.concatenate((y_eq[:, 4], y_input1[:, 4], y_input2[:, 4]))
    Dimer = AB_values + ApB_values

    return t, A_values, B_values, Ap_values, AB_values, ApB_values, eq_duration, input1_duration, input2_duration, inputs


def format_ode_plot(inputs, eq_duration, input1_duration, input2_duration, modeltype, legend):
    if inputs == 2:
        pass
    elif inputs == 1:
        plt.xlim(0, eq_duration+input1_duration)
    else:
        pass

    if modeltype == "all_concentrations":
        plt.ylabel("Concentration (M)", size=15)
    else:
        plt.ylabel("Dimer fraction", size=15)
        plt.ylim(0, 1)
    plt.margins(0)
    plt.grid(axis="y")
    plt.xlabel("Time", size=15)
    plt.axvline(x=eq_duration, color="grey", alpha=0.5, zorder=-99)
    plt.axvline(x=eq_duration+input1_duration,
                color="grey", alpha=0.5, zorder=-100)
    plt.axvspan(0, eq_duration, facecolor="grey", alpha=0.1, zorder=-101)
    plt.axvspan(eq_duration, eq_duration+input1_duration,
                facecolor="red", alpha=0.1, zorder=-102)
    plt.axvspan(eq_duration+input1_duration, eq_duration+input1_duration+input2_duration,
                facecolor="blue", alpha=0.1, zorder=-103)

    plt.xlabel("Time (seconds)")
    if legend == True:
        leg = plt.legend()
        leg.set_draggable(state=True)
    else:
        pass
    plt.show()


def ode_switch_model_plot(y0: dict, kinetic_params: dict, eq_duration=200, input1_duration=600, input2_duration=600, inputs=2, pp_background="without_pp_background", modeltype="dimer_fraction", figsize=(6.4, 4.8), legend: bool = True):
    data = ode_switch_model(y0=y0, kinetic_params=kinetic_params, eq_duration=eq_duration,
                            input1_duration=input1_duration, inputs=inputs, pp_background=pp_background)
    if not modeltype in [None, "all_concentrations", "dimer_fraction"]:
        raise ValueError(
            "Graph type must be 'all_concentrations', 'dimer_fraction'.")

    dimer = data[4] + data[5]
    dimer_fraction = dimer/min(y0["A"], y0["B"])

    plt.figure(figsize=figsize)

    t, A_values, B_values, Ap_values, AB_values, ApB_values = data[
        0], data[1], data[2], data[3], data[4], data[5]
    if modeltype == "all_concentrations":
        plt.plot(t, A_values, label="A", color="Blue", linewidth=1)
        plt.plot(t, Ap_values, label="Ap", color="Blue",
                 linewidth=1, linestyle="--")
        plt.plot(t, B_values, label="B", color="orange", linewidth=1)
        plt.plot(t, AB_values, label="AB", color="red", linewidth=1)
        plt.plot(t, ApB_values, label="ApB", color="red",
                 linewidth=0.8, linestyle="--")
        plt.ylabel("Concentration (M)", size=15)
        plt.ticklabel_format(useOffset=False, useMathText=True)

    elif modeltype == "dimer_fraction":
        plt.plot(t, dimer_fraction, label="Dimer fraction",
                 color="black", linewidth=2.5)

    format_ode_plot(inputs, eq_duration, input1_duration,
                    input2_duration, modeltype, legend)


def param_sweep_plot(y0: dict, kinetic_params: dict, param_values: list, variable_param: str, eq_duration=200, input1_duration=600, input2_duration=600, pp_background="without_pp_background", inputs=2, figsize=(6.4, 4.8), legend: bool = True):
    if variable_param not in ["kon1", "koff1", "kon2", "koff2", "cat_kin1", "cat_pp1", "cat_kin2", "cat_pp2"]:
        raise ValueError(
            "variable_param must be a valid parameter: 'kon1', 'koff1', 'kon2', 'koff2', 'cat_kin1', 'cat_pp1', 'cat_kin2', 'cat_pp2'")

    # Normalize the log10 of parameter values
    norm_param_values = [(np.log10(val) - np.log10(min(param_values))) /
                         (np.log10(max(param_values)) - np.log10(min(param_values))) for val in param_values]

    plt.figure(figsize=figsize)

    for i, value in enumerate(param_values):
        kinetic_params[variable_param] = value
        data = ode_switch_model(y0=y0, kinetic_params=kinetic_params, eq_duration=eq_duration, input1_duration=input1_duration,
                                input2_duration=input2_duration, pp_background=pp_background, inputs=2)
        dimer = data[4] + data[5]
        dimer_fraction = dimer/min(y0["A"], y0["B"])
        # Use the colormap to color the line based on the normalized log10 of the parameter value
        color = cm.viridis(norm_param_values[i])
        plt.plot(data[0], dimer_fraction,
                 label=f'{variable_param}={value:.2e}', color=color)

    format_ode_plot(inputs, eq_duration, input1_duration,
                    input2_duration, modeltype="dimer_fraction", legend=legend)


def multi_param_sweep_plot(y0: dict, kinetic_params: dict, param_values: list, variable_params: tuple, eq_duration=200, input1_duration=600, input2_duration=600, pp_background="without_pp_background", inputs=2, figsize=(6.4, 4.8), legend: bool = True):
    var1, var2 = variable_params
    if var1 not in ["kon1", "koff1", "kon2", "koff2", "cat_kin1", "cat_pp1", "cat_kin2", "cat_pp2"] or var2 not in ["kon1", "koff1", "kon2", "koff2", "cat_kin1", "cat_pp1", "cat_kin2", "cat_pp2"]:
        raise ValueError(
            "variable_params must be a tuple of two parameters, one of which must be one of the following: 'kon1', 'koff1', 'kon2', 'koff2', 'cat_kin1', 'cat_pp1', 'cat_kin2', 'cat_pp2'")

    # Normalize the log10 of parameter values
    norm_param_values = [(np.log10(val[0]) - np.log10(min([i[0] for i in param_values]))) /
                         (np.log10(max([i[0] for i in param_values])) - np.log10(min([i[0] for i in param_values]))) for val in param_values]

    plt.figure(figsize=figsize)

    for i, value in enumerate(param_values):
        kinetic_params[var1] = value[0]
        kinetic_params[var2] = value[1]
        data = ode_switch_model(y0=y0, kinetic_params=kinetic_params, eq_duration=eq_duration, input1_duration=input1_duration,
                                input2_duration=input2_duration, pp_background=pp_background, inputs=2)
        dimer = data[4] + data[5]
        dimer_fraction = dimer/min(y0["A"], y0["B"])
        # Use the colormap to color the line based on the normalized log10 of the first parameter value
        color = cm.viridis(norm_param_values[i])
        plt.plot(data[0], dimer_fraction,
                 label=f'{var1}={value[0]:.2e}, {var2}={value[1]:.2e}', color=color)

    format_ode_plot(inputs, eq_duration, input1_duration,
                    input2_duration, modeltype="dimer_fraction", legend=legend)


def kon_koff_sweep_plot(y0: dict, kinetic_params: dict, param_tuples1: list, param_tuples2: list, eq_duration=200, input1_duration=600, input2_duration=600, pp_background="without_pp_background", inputs=2, figsize=(6.4, 4.8), legend: bool = True):
    if len(param_tuples1) != len(param_tuples2):
        raise ValueError("Both lists must have the same length.")

    plt.figure(figsize=figsize)

    for i, (tuple1, tuple2) in enumerate(zip(param_tuples1, param_tuples2)):
        kinetic_params["kon1"] = tuple1[0]
        kinetic_params["koff1"] = tuple1[1]
        kinetic_params["kon2"] = tuple2[0]
        kinetic_params["koff2"] = tuple2[1]
        data = ode_switch_model(y0=y0, kinetic_params=kinetic_params, eq_duration=eq_duration, input1_duration=input1_duration,
                                input2_duration=input2_duration, pp_background=pp_background, inputs=2)
        dimer = data[4] + data[5]
        dimer_fraction = dimer/min(y0["A"], y0["B"])
        color = cm.viridis(i/len(param_tuples1))
        plt.plot(data[0], dimer_fraction,
                 label=f'kon1={tuple1[0]:.2e}, koff1={tuple1[1]:.2e}, kon2={tuple2[0]:.2e}, koff2={tuple2[1]:.2e}', color=color)

    format_ode_plot(inputs, eq_duration, input1_duration,
                    input2_duration, modeltype="dimer_fraction", legend=legend)


def generate_equal_tuples(min_val: float, max_val: float, num_values: int):
    log_min_val = np.log10(min_val)
    log_max_val = np.log10(max_val)
    log_vals = np.logspace(log_min_val, log_max_val, num_values)
    param_tuples = [(val, val) for val in log_vals]
    return param_tuples


def generate_kd_tuples(ratio: float, min_val: float, max_val: float, num_values: int):
    log_min_val = np.log10(min_val)
    log_max_val = np.log10(max_val)
    log_vals = np.logspace(log_min_val, log_max_val, num_values)
    param_tuples = [(val, val*ratio) for val in log_vals]
    return param_tuples


y0 = {"A": 1e-6, "B": 1e-6, "Ap": 0, "AB": 0, "ApB": 0}
kinetic_params = {"kon1": 1e4, "koff1": 2.25e-2, "kon2": 1e4, "koff2": 8.30e-4,
                  "cat_kin1": 1e-2, "cat_pp1": 0, "cat_kin2": 1e-2, "cat_pp2": 0}
variable_params = ("cat_kin2", "cat_pp2")
param_values = generate_equal_tuples(1e-4, 1e-2, 7)

multi_param_sweep_plot(y0, kinetic_params=kinetic_params, variable_params=variable_params,
                       param_values=param_values, inputs=2, input1_duration=900, input2_duration=800, figsize=(9.6, 4.8))

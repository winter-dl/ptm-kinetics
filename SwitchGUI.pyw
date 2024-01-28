import PySimpleGUI as sg
import eqmodel as eq

sg.theme("SystemDefault")

layout=[

[sg.Button("Show me an example!", size=(20,1), key="EXAMPLE")],

[sg.Text(f"Enter the dissociation constants (in M) between the two subunits of the protein switch:")],

[
    sg.Text("Kd value in the ABSENCE of input: Kd1 ="),
    sg.Input(size=(8, 1), key="KD_WO_INPUT"),
    sg.Text("M")
    ],

[
    sg.Text("Kd value in the PRESENCE of input: Kd2 ="),
    sg.Input(size=(8, 1), enable_events=True, key="KD_WITH_INPUT"),
    sg.Text("M")
    ],

[sg.HSeparator(pad=(0, 10))],

[
    sg.Text(f"Below, you can enter more than one value for protein switch concentrations (in M) or dimer/monomer baselines. Separate the values by commas.")
    ],

[sg.Text("Let's switch!")],

[sg.HSeparator(pad=(0, 10))],

[
    sg.Text("Concentration(s):"),
    sg.Input(size=(40, 1), enable_events=True, key="CONCENTRATIONS"),
    sg.Text(f"M")
    ],

[sg.Button("Switch behaviour at chosen concentration(s)", size=(24,2), pad=(0,5), key="SWITCH"), sg.Text("The overall behaviour is calculated from Kd1 and Kd2. Details are plotted for the chosen concentration(s).")],

[sg.HSeparator(pad=(0, 10))],

[
    sg.Text("Baseline value(s):"),
    sg.Input(size=(40, 1), enable_events=True, key="BASELINES"),
    sg.Text(f"(e.g., 0.03)")
    ],

[sg.Button(f"Switch behaviour at chosen\nDIMER baselines(s)", size=(24,2), pad=(0,5), key="DIMERBASELINE"), sg.Text(f"The overall behaviour is calculated from Kd1 and Kd2. Details are plotted for the concentration(s) required\nto produce the chosen DIMER fraction baseline(s).")],

[sg.Button(f"Switch behaviour at chosen\nMONOMER baselines(s)", size=(24,2), pad=(0,5), key="MONOMERBASELINE"), sg.Text(f"The overall behaviour is calculated from Kd1 and Kd2. Details are plotted for the concentration(s) required\nto produce the chosen MONOMER fraction baseline(s).")],

[sg.HSeparator(pad=(0, 10))],

[sg.Button("Optimise for dynamic range", size=(24,2), pad=(0,5), key="C_OPT"), sg.Text(f"The overall behaviour is calculated from Kd1 and Kd2. Details are plotted for the concentration at which\nthe physicochemical input causes the greatest absolute change in dimerisation.")],

[sg.Checkbox(f"Show annotations.", default=True, key="ANNO"), sg.Text(f"Annotations are draggable. Unchecking this box yields clean graphs where simple annotations can be brought up\nby clicking on points.")],

[sg.Output(size=(100, 15), background_color="black", text_color='white', font=("Courier", 10))],

]

window = sg.Window("Dimeric Switch", layout)

while True:
    event, value = window.read() 
    if event in (sg.WIN_CLOSED, 'Exit'):
        break
    elif event == "EXAMPLE":
        window["KD_WO_INPUT"]("300e-9")
        window["KD_WITH_INPUT"]("20e-9")
        window["CONCENTRATIONS"]("10e-9")
        window["BASELINES"]("0.01, 0.05")
    elif event == "SWITCH":
        try:
            kd1 = float(value["KD_WO_INPUT"])
            kd2 = float(value["KD_WITH_INPUT"])
            clist = value["CONCENTRATIONS"].split(",")
            c=[]
            for i in clist:
                c.append(float(i))
            switch = eq.DimericSwitch(kd1, kd2)  
            if value["ANNO"] == True:
                switch.get_behaviour(*c, annotations=True)
            else:
                switch.get_behaviour(*c, annotations=False)
        except:
            sg.Popup("Kd1, Kd2, and concentration values must be numbers.", title="Error")
    elif event == "DIMERBASELINE":
        try:
            kd1 = float(value["KD_WO_INPUT"])
            kd2 = float(value["KD_WITH_INPUT"])
            blist = value["BASELINES"].split(",")
            b=[]
            for i in blist:
                b.append(float(i))
            switch = eq.DimericSwitch(kd1, kd2)
            for i in b:
                switch.get_key_values_at_dimer_baseline(i)
                print(f"")
            if value["ANNO"] == True:
                switch.get_behaviour(*b, type='dimerbaseline', annotations=True)
            else:
                switch.get_behaviour(*b, type='dimerbaseline', annotations=False)
            print(f'{100*"-"}\n')
        except:
            sg.Popup("Kd1, Kd2, and baseline values must be numbers.", title="Error")
    elif event == "MONOMERBASELINE":
        try:
            kd1 = float(value["KD_WO_INPUT"])
            kd2 = float(value["KD_WITH_INPUT"])
            blist = value["BASELINES"].split(",")
            b=[]
            for i in blist:
                b.append(float(i))
            switch = eq.DimericSwitch(kd1, kd2)
            for i in b:
                switch.get_key_values_at_monomer_baseline(i)
                print(f"")
            if value["ANNO"] == True:
                switch.get_behaviour(*b, type='monomerbaseline', annotations=True)
            else:
                switch.get_behaviour(*b, type='monomerbaseline', annotations=False)
            print(f'{100*"-"}\n')
        except:
            sg.Popup("Kd1, Kd2, and baseline values must be numbers.", title="Error")
    elif event == "C_OPT":
        try:
            kd1 = float(value["KD_WO_INPUT"])
            kd2 = float(value["KD_WITH_INPUT"])
            c = eq.optimise_concentration(kd1, kd2)["Optimal concentration"]
            window["CONCENTRATIONS"]("{:.3e}".format(c))
            switch = eq.DimericSwitch(kd1, kd2)
            switch.get_key_values()
            if value["ANNO"] == True:
                switch.get_behaviour(type="optimal", annotations=True)
            else:
                switch.get_behaviour(type="optimal", annotations=False)
            print(f'\n{100*"-"}\n')
        except:
            sg.Popup("Kd1 and Kd2 values must be numbers.", title="Error")
window.close()
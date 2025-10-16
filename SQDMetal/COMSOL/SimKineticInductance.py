from SQDMetal.COMSOL.Model import COMSOL_Simulation_Base

import jpype.types as jtypes
import geopandas as gpd
import shapely
import matplotlib.pyplot as plt
import numpy as np

class COMSOL_Simulation_KinInd(COMSOL_Simulation_Base):

    def __init__ (self, model, sim_type="eigfreq"):
        self.sqdmodel = model
        if sim_type not in ["eigfreq","sparams"]:
            raise ValueError("sim_type must be 'eigfreq' or 'sparams'")
        self.sim_type = sim_type
        self.model = model._get_java_comp()

    ## Customized part
    def import_parameters(self, param_file="mph/Al_params.txt"):
        self.model.param().loadFile(param_file)

    def redefine_geometry(self, chip_dim_x="0.003", chip_dim_y="0.003"):
        self.model.component("comp1").geom("geom1").feature("blk_chip").set("size", jtypes.JArray(jtypes.JString)([chip_dim_x, chip_dim_y, "t_sub"]))
        self.model.component("comp1").geom("geom1").feature("blk_chip").set("pos", jtypes.JArray(jtypes.JString)(["0", "0", "-t_sub/2"]))
        self.model.component("comp1").geom("geom1").feature("blk_boundary").set("size", jtypes.JArray(jtypes.JString)([chip_dim_x, chip_dim_y, "t_air"]))
        self.model.component("comp1").geom("geom1").feature("blk_boundary").set("pos", jtypes.JArray(jtypes.JString)(["0", "0", "t_air/2"]))
        self.model.component("comp1").geom("geom1").runPre("fin")
        self.model.component("comp1").geom("geom1").run("selFused1")
        self.model.component("comp1").geom("geom1").feature().create("ext1", "Extrude")
        self.model.component("comp1").geom("geom1").feature("ext1").set("workplane", "wp1")
        self.model.component("comp1").geom("geom1").feature("ext1").selection("input").set("wp1")
        self.model.component("comp1").geom("geom1").feature("ext1").setIndex("distance", "t", 0)
        self.model.component("comp1").geom("geom1").runPre("fin")
        self.model.component("comp1").geom("geom1").run()
        self.model.component("comp1").geom("geom1").runPre("ext1")
        self.model.component("comp1").geom("geom1").run("selFused1")
        self.model.component("comp1").geom("geom1").create("difsel1", "DifferenceSelection")
        self.model.component("comp1").geom("geom1").feature("difsel1").set("entitydim", jtypes.JInt(2))
        self.model.component("comp1").geom("geom1").feature("difsel1").set("add", jtypes.JArray(jtypes.JString)(["wp1_selcond0"]))
        self.model.component("comp1").geom("geom1").feature("difsel1").set("subtract", jtypes.JArray(jtypes.JString)(["selFused0", "selFused1"]))
        self.model.component("comp1").geom("geom1").runPre("fin")
        self.model.component("comp1").geom("geom1").feature("difsel1").set("entitydim", jtypes.JInt(3))
        self.model.component("comp1").geom("geom1").runPre("fin")
        self.model.component("comp1").geom("geom1").feature("difsel1").set("entitydim", jtypes.JInt(2))
        self.model.component("comp1").geom("geom1").runPre("fin")
        self.model.component("comp1").geom("geom1").runPre("difsel1")
        self.model.component("comp1").geom("geom1").run("selFused1")
        self.model.component("comp1").geom("geom1").create("unisel1", "UnionSelection")
        self.model.component("comp1").geom("geom1").feature("unisel1").label("condAll")
        self.model.component("comp1").geom("geom1").feature("unisel1").set("entitydim", jtypes.JInt(2))
        self.model.component("comp1").geom("geom1").feature("unisel1").set("input", jtypes.JArray(jtypes.JString)(["selFused0", "selFused1"]))
        self.model.component("comp1").geom("geom1").runPre("difsel1")
        self.model.component("comp1").geom("geom1").feature("wp1").geom().run("difGND")
        self.model.component("comp1").geom("geom1").feature("wp1").geom().create("r1", "Rectangle")
        self.model.component("comp1").geom("geom1").feature("wp1").geom().feature("r1").set("size", jtypes.JArray(jtypes.JString)([chip_dim_x, chip_dim_y]))
        self.model.component("comp1").geom("geom1").feature("wp1").geom().feature("r1").set("base", "center")
        self.model.component("comp1").geom("geom1").feature("wp1").geom().feature("r1").set("pos", jtypes.JArray(jtypes.JInt)([0, 0]))
        self.model.component("comp1").geom("geom1").feature("wp1").geom().run("r1")
        self.model.component("comp1").geom("geom1").feature("wp1").geom().create("sel1", "ExplicitSelection")
        self.model.component("comp1").geom("geom1").feature("wp1").geom().feature("sel1").selection("selection").set("r1", 1)
        self.model.component("comp1").geom("geom1").runPre("difsel1")
        self.model.component("comp1").geom("geom1").feature("difsel1").set("add", jtypes.JArray(jtypes.JString)(["wp1_sel1"]))
        self.model.component("comp1").geom("geom1").feature("difsel1").set("subtract", jtypes.JArray(jtypes.JString)(["unisel1"]))
        self.model.component("comp1").geom("geom1").runPre("fin")
        self.model.component("comp1").geom("geom1").run()

    def set_materials(self, sub=[1], metal=[2,4]):
        self.sqdmodel._create_material('Vacuum', 1.0, 1.0, 0.0)
        self.model.component("comp1").material('Vacuum').selection().all()
        self.sqdmodel._create_material('Substrate', 0.0, 1.0, 0.0)
        self.sqdmodel._create_material('Metal', 0.0, 1.0, 3.77e7)
        self.model.component("comp1").material("Substrate").selection().set(jtypes.JArray(jtypes.JInt)(sub))
        self.model.component("comp1").material("Metal").selection().set(jtypes.JArray(jtypes.JInt)(metal))
        self.model.component("comp1").material("Substrate").propertyGroup("def").set("relpermittivity", jtypes.JArray(jtypes.JString)(["eps_Si"]))
        self.model.component("comp1").material("Metal").propertyGroup("def").set("electricconductivity", jtypes.JArray(jtypes.JString)(["-1j*alpha/emw.freq"]))
        self.model.component("comp1").physics().create("emw", "ElectromagneticWaves", "geom1")

    def mesh(self):
        if self.sim_type == "eigfreq":
            self.model.component("comp1").mesh("mesh1").create("ftri1", "FreeTri")
            self.model.component("comp1").mesh("mesh1").feature("ftri1").selection().set(9, 17, 22)
            self.model.component("comp1").mesh("mesh1").feature("ftri1").create("size1", "Size")
            self.model.component("comp1").mesh("mesh1").feature("ftri1").feature("size1").set("hauto",jtypes.JInt(1))
            self.model.component("comp1").mesh("mesh1").feature("ftri1").feature("size1").set("custom", True)
            self.model.component("comp1").mesh("mesh1").feature("ftri1").feature("size1").set("hmaxactive", True)
            self.model.component("comp1").mesh("mesh1").feature("ftri1").feature("size1").set("hmax", 6.0E-4)
            self.model.component("comp1").mesh("mesh1").feature("ftri1").feature("size1").set("hminactive", True)
            self.model.component("comp1").mesh("mesh1").feature("ftri1").feature("size1").set("hcurveactive", False)
            self.model.component("comp1").mesh("mesh1").feature("ftri1").feature("size1").set("hnarrowactive", True)
            self.model.component("comp1").mesh("mesh1").create("swe1", "Sweep")
            self.model.component("comp1").mesh("mesh1").feature("swe1").selection().geom("geom1", jtypes.JInt(3))
            self.model.component("comp1").mesh("mesh1").feature("swe1").selection().set(2, 4, 5)
            self.model.component("comp1").mesh("mesh1").feature("swe1").selection("sourceface").set(9, 17, 22)
            self.model.component("comp1").mesh("mesh1").feature("swe1").selection("targetface").set(6, 16, 21)
            self.model.component("comp1").mesh("mesh1").feature("swe1").create("dis1", "Distribution")
            self.model.component("comp1").mesh("mesh1").feature("swe1").feature("dis1").set("numelem", jtypes.JInt(1))
            self.model.component("comp1").mesh("mesh1").create("ftet1", "FreeTet")
            self.model.component("comp1").mesh("mesh1").feature("ftet1").selection().geom("geom1", jtypes.JInt(3))
            self.model.component("comp1").mesh("mesh1").feature("ftet1").selection().set(jtypes.JInt(3))
            self.model.component("comp1").mesh("mesh1").feature().duplicate("ftet2", "ftet1")
            self.model.component("comp1").mesh("mesh1").feature("ftet2").selection().set(jtypes.JInt(1))
            self.model.component("comp1").mesh("mesh1").feature("ftet2").selection().remaining()
        else:
            self.model.component("comp1").mesh("mesh1").create("ftri1", "FreeTri")
            self.model.component("comp1").mesh("mesh1").feature("ftri1").selection().set(9, 14, 18, 22, 26, 36)
            self.model.component("comp1").mesh("mesh1").create("swe1", "Sweep")
            self.model.component("comp1").mesh("mesh1").feature("swe1").selection("sourceface").set(9, 14, 18, 22, 26, 36)
            self.model.component("comp1").mesh("mesh1").feature("swe1").selection().geom("geom1", 3)
            self.model.component("comp1").mesh("mesh1").feature("swe1").selection().set(2, 4, 5, 6, 7, 8)
            self.model.component("comp1").mesh("mesh1").feature("swe1").selection("targetface").set(6, 13, 17, 21, 25, 35)
            self.model.component("comp1").mesh("mesh1").feature("ftri1").create("size1", "Size")
            self.model.component("comp1").mesh("mesh1").feature("ftri1").feature("size1").set("hauto", jtypes.JInt(1))
            self.model.component("comp1").mesh("mesh1").feature("ftri1").feature("size1").set("custom", True)
            self.model.component("comp1").mesh("mesh1").feature("ftri1").feature("size1").set("hmaxactive", True)
            self.model.component("comp1").mesh("mesh1").feature("ftri1").feature("size1").set("hminactive", True)
            self.model.component("comp1").mesh("mesh1").feature("ftri1").feature("size1").set("hmax", 6.0E-4)
            self.model.component("comp1").mesh("mesh1").feature("ftri1").feature("size1").set("hcurveactive", True)
            self.model.component("comp1").mesh("mesh1").feature("ftri1").feature("size1").set("hnarrowactive", True)
            self.model.component("comp1").mesh("mesh1").feature("swe1").create("dis1", "Distribution")
            self.model.component("comp1").mesh("mesh1").feature("swe1").feature("dis1").set("numelem", jtypes.JInt(1))
            self.model.component("comp1").mesh("mesh1").create("ftet1", "FreeTet")
            self.model.component("comp1").mesh("mesh1").feature("ftet1").selection().geom("geom1", 3)
            self.model.component("comp1").mesh("mesh1").feature("ftet1").selection().set(3)
            self.model.component("comp1").mesh("mesh1").create("ftet2", "FreeTet")

            self.model.component("comp1").physics("emw").create("lport1", "LumpedPort", 2)
            self.model.component("comp1").physics("emw").feature("lport1").selection().set(15)
            self.model.component("comp1").physics("emw").create("lport2", "LumpedPort", 2)
            self.model.component("comp1").physics("emw").feature("lport2").selection().set(584)
            self.model.component("comp1").physics("emw").feature("lport2").set("PortExcitation", "off")
        
        self.model.component("comp1").mesh("mesh1").run()

    def setup_study(self, search_around_freq_GHz="3", linearization_freq_GHz="4", s_params=False):
        self.model.study().create("std1")
        self.model.study("std1").create("eig", "Eigenfrequency")
        self.model.study("std1").feature("eig").set("conrad", "1")
        self.model.study("std1").feature("eig").set("linpsolnum", "auto")
        self.model.study("std1").feature("eig").set("solnum", "auto")
        self.model.study("std1").feature("eig").set("notsolnum", "auto")
        self.model.study("std1").feature("eig").set("ngenAUX", "1")
        self.model.study("std1").feature("eig").set("goalngenAUX", "1")
        self.model.study("std1").feature("eig").set("ngenAUX", "1")
        self.model.study("std1").feature("eig").set("goalngenAUX", "1")
        self.model.study("std1").feature("eig").setSolveFor("/physics/emw", True)
        self.model.study("std1").feature("eig").set("neigsactive", True)
        self.model.study("std1").feature("eig").set("neigs", jtypes.JInt(1))
        self.model.study("std1").feature("eig").set("eigwhich", "lr")
        self.model.sol().create("sol1")
        self.model.sol("sol1").study("std1")
        self.model.study("std1").feature("eig").set("notlistsolnum", jtypes.JInt(1))
        self.model.study("std1").feature("eig").set("notsolnum", "auto")
        self.model.study("std1").feature("eig").set("listsolnum", jtypes.JInt(1))
        self.model.study("std1").feature("eig").set("solnum", "auto")
        self.model.study("std1").feature("eig").set("linplistsolnum", jtypes.JArray(jtypes.JString)(["1"]))
        self.model.study("std1").feature("eig").set("linpsolnum", "auto")
        self.model.study("std1").feature("eig").set("shift", jtypes.JString(f"{search_around_freq_GHz}[GHz]"))
        self.model.sol("sol1").create("st1", "StudyStep")
        self.model.sol("sol1").feature("st1").set("study", "std1")
        self.model.sol("sol1").feature("st1").set("studystep", "eig")
        self.model.sol("sol1").create("v1", "Variables")
        self.model.sol("sol1").feature("v1").set("control", "eig")
        self.model.sol("sol1").create("e1", "Eigenvalue")
        self.model.sol("sol1").feature("e1").set("control", "eig")
        self.model.sol("sol1").feature("e1").feature("aDef").set("complexfun", True)
        self.model.sol("sol1").feature("e1").create("d1", "Direct")
        self.model.sol("sol1").feature("e1").feature("d1").set("linsolver", "pardiso")
        self.model.sol("sol1").feature("e1").feature("d1").label("Suggested Direct Solver (emw)")
        self.model.sol("sol1").feature("e1").create("i1", "Iterative")
        self.model.sol("sol1").feature("e1").feature("i1").set("linsolver", "gmres")
        self.model.sol("sol1").feature("e1").feature("i1").set("prefuntype", "right")
        self.model.sol("sol1").feature("e1").feature("i1").set("itrestart", "300")
        self.model.sol("sol1").feature("e1").feature("i1").label("Suggested Iterative Solver (emw)")
        self.model.sol("sol1").feature("e1").feature("i1").create("mg1", "Multigrid")
        self.model.sol("sol1").feature("e1").feature("i1").feature("mg1").set("iter", "1")
        self.model.sol("sol1").feature("e1").feature("i1").feature("mg1").feature("pr").create("sv1", "SORVector")
        self.model.sol("sol1").feature("e1").feature("i1").feature("mg1").feature("pr").feature("sv1").set("prefun", "sorvec")
        self.model.sol("sol1").feature("e1").feature("i1").feature("mg1").feature("pr").feature("sv1").set("iter", jtypes.JInt(2))
        self.model.sol("sol1").feature("e1").feature("i1").feature("mg1").feature("pr").feature("sv1").set("relax", jtypes.JInt(1))
        self.model.sol("sol1").feature("e1").feature("i1").feature("mg1").feature("pr").feature("sv1").set("sorvecdof", jtypes.JArray(jtypes.JString)(["comp1_E"]))
        self.model.sol("sol1").feature("e1").feature("i1").feature("mg1").feature("po").create("sv1", "SORVector")
        self.model.sol("sol1").feature("e1").feature("i1").feature("mg1").feature("po").feature("sv1").set("prefun", "soruvec")
        self.model.sol("sol1").feature("e1").feature("i1").feature("mg1").feature("po").feature("sv1").set("iter", jtypes.JInt(2))
        self.model.sol("sol1").feature("e1").feature("i1").feature("mg1").feature("po").feature("sv1").set("relax", jtypes.JInt(1))
        self.model.sol("sol1").feature("e1").feature("i1").feature("mg1").feature("po").feature("sv1").set("sorvecdof", jtypes.JArray(jtypes.JString)(["comp1_E"]))
        self.model.sol("sol1").feature("e1").feature("i1").feature("mg1").feature("cs").create("d1", "Direct")
        self.model.sol("sol1").feature("e1").feature("i1").feature("mg1").feature("cs").feature("d1").set("linsolver", "pardiso")
        self.model.sol("sol1").feature("e1").feature("d1").active(True)
        self.model.sol("sol1").attach("std1")
        self.model.sol("sol1").feature("e1").set("eigref", linearization_freq_GHz)
        self.model.sol("sol1").feature("e1").feature("d1").set("linsolver", "mumps")

        if s_params:
            pass

    def run(self, s_params=False):

        # --- Creazione risultati (come nel tuo codice) ---
        self.model.result().create("pg1", "PlotGroup3D")
        self.model.result("pg1").label("Electric Field (emw)")
        self.model.result("pg1").set("frametype", "spatial")
        self.model.result("pg1").set("showlegendsmaxmin", True)
        self.model.result("pg1").set("defaultPlotID", "ElectromagneticWaves/phys1/pdef1/pcond1/pg1")
        self.model.result("pg1").feature().create("mslc1", "Multislice")
        self.model.result("pg1").feature("mslc1").label("Multislice")
        self.model.result("pg1").feature("mslc1").set("smooth", "internal")
        self.model.result("pg1").feature("mslc1").set("data", "parent")
        self.model.result("pg1").feature("mslc1").feature().create("filt1", "Filter")
        self.model.result("pg1").feature("mslc1").feature("filt1").set("expr", "!isScalingSystemDomain")

        self.model.result().numerical().create("gev1", "EvalGlobal")
        self.model.result().numerical("gev1").label("Eigenfrequencies (emw)")
        self.model.result().numerical("gev1").set("data", "dset1")
        self.model.result().numerical("gev1").set("expr", jtypes.JArray(jtypes.JString)(["emw.freq", "emw.Qfactor"]))
        self.model.result().numerical("gev1").set("unit", jtypes.JArray(jtypes.JString)(["GHz", "1"]))
        self.model.result().table().create("tbl1", "Table")
        self.model.result().numerical("gev1").set("table", "tbl1")

        # --- Esecuzione solver ---
        self.model.sol("sol1").runAll()

        # --- Calcolo risultati numerici ---
        self.model.result().numerical("gev1").setResult()

        # --- Estrazione dati dalla tabella ---
        data = self.model.result().numerical("gev1").getReal()  # oppure getData(), dipende dalla versione COMSOL

        # 'data' è una lista di liste [[freq1, Q1], [freq2, Q2], ...]
        first_freq = None
        if data and len(data) > 0:
            first_freq = float(data[0][0])  # prima autofrequenza (in GHz)
            print(f"Prima autofrequenza trovata: {first_freq} GHz")
        else:
            print("⚠️ Nessuna autofrequenza trovata!")

        if s_params:
            pass

        return first_freq

    def plot_results(self):
        self.model.result("pg1").run()
        self.model.result("pg1").feature("mslc1").set("multiplanexmethod", "number")
        self.model.result("pg1").feature("mslc1").set("xnumber", "0")
        self.model.result("pg1").feature("mslc1").set("ynumber", "0")
        self.model.result("pg1").feature("mslc1").set("multiplanezmethod", "coord")
        self.model.result("pg1").feature("mslc1").set("zcoord", jtypes.JInt(0))
        self.model.result("pg1").feature("mslc1").set("colortable", "ThermalWaveDark")
        self.model.result("pg1").feature("mslc1").create("def1", "Deform")
        self.model.result("pg1").run()
        self.model.result("pg1").feature("mslc1").feature("def1").set("expr", jtypes.JArray(jtypes.JString)(["0", "0", "emw.normE"]))
        self.model.result("pg1").run()
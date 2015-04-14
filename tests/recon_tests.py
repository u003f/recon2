"""
Various tests to compare different flavours of recon 2.
Updated: 10 Apr 15 by Kieran Smallbone
"""

import numpy
import os
import scipy.sparse as sparse

# http://www.gurobi.com/documentation/6.0/quickstart_mac/py_python_interface.html
import gurobipy
# http://frank-fbergmann.blogspot.co.uk/2014/05/libsbml-python-bindings-5101.html
import libsbml


INF = numpy.inf
NAN = numpy.nan


def run_all():
    """Runs all the tests on all the models"""

    model_names, model_path = list_models()
    model_names = ['recon_2.2']
    for name in model_names:
        print '\n%s'%name
        filename = os.path.join(model_path, name + '.xml')
        max_fluxes(filename)


def list_models():
    """
    model_names, model_path = list_models()
    returns
    model_names: list of SBML models in the directory ../models
    model_path: the full path to the directory ../models
    """

    tests_path = os.path.dirname(__file__)
    model_path = os.path.join(tests_path, '..', 'models')
    model_path = os.path.normpath(model_path)
    model_names = []
    for filename in os.listdir(model_path):
        shortname, extension = os.path.splitext(filename)
        if extension == '.xml':
            model_names.append(shortname)
    model_names.sort(reverse=True)  # ~ most recent first
    return model_names, model_path


def max_fluxes(model_filename):
    """
    Written to mimic neilswainston matlab function maxFluxes
    """

    media = [
        'EX_ca2(e)',
        'EX_cl(e)',
        'EX_fe2(e)',
        'EX_fe3(e)',
        'EX_h(e)',
        'EX_h2o(e)',
        'EX_k(e)',
        'EX_na1(e)',
        'EX_nh4(e)',
        'EX_so4(e)',
        'EX_pi(e)'
        ]
    objective = 'DM_atp_c_'
    for normoxic in [True, False]:
        for carbon_source in [
                # sugars
                'EX_glc(e)',
                'EX_fru(e)',
                # fatty acids
                'EX_ppa(e)',        # C3:0
                'EX_but(e)',        # C4:0
                'EX_hx(e)',         # C6:0
                'EX_octa(e)',       # C8:0
                'EX_HC02175(e)',    # C10:0
                'EX_HC02176(e)',    # C12:0
                'EX_ttdca(e)',      # C14:0
                'EX_hdca(e)',       # C16:0
                'EX_ocdca(e)',      # C18:0
                'EX_arach(e)',      # C20:0
                'EX_docosac_',      # C22:0
                'EX_lgnc(e)',       # C24:0
                'EX_hexc(e)',       # C26:0
                # amino acids
                'EX_ala_L(e)',
                'EX_arg_L(e)',
                'EX_asn_L(e)',
                'EX_asp_L(e)',
                'EX_cys_L(e)',
                'EX_gln_L(e)',
                'EX_glu_L(e)',
                'EX_gly(e)',
                'EX_his_L(e)',
                'EX_ile_L(e)',
                'EX_leu_L(e)',
                'EX_lys_L(e)',
                'EX_met_L(e)',
                'EX_phe_L(e)',
                'EX_pro_L(e)',
                'EX_ser_L(e)',
                'EX_thr_L(e)',
                'EX_trp_L(e)',
                'EX_tyr_L(e)',
                'EX_val_L(e)',
                ]:
            f_opt = max_flux(model_filename, carbon_source, objective, normoxic, media)
            print '%s:\t%g'%(carbon_source, f_opt)


def max_flux(model_filename, carbon_source, objective, normoxic, media):
    """
    Written to mimic neilswainston matlab function maxFlux
    """
    sbml = read_sbml(model_filename)
    set_infinite_bounds(sbml)
    # block import reactions
    block_all_imports(sbml)
    # define carbon source
    set_import_bounds(sbml, carbon_source, 1)
    # define media
    set_import_bounds(sbml, media, INF)
    if normoxic:
        set_import_bounds(sbml, 'EX_o2(e)', INF)
    # specify objective and maximise
    change_objective(sbml, objective)
    # avoid infinities
    obj_max = 1e6
    change_rxn_bounds(sbml, objective, obj_max, 'u')
    v_sol, f_opt = optimize_cobra_model(sbml)
    if f_opt > 0.9 * obj_max:
        f_opt = INF
    return f_opt


def read_sbml(filename):
    """
    Read an SBML file from specified path.
    Copied from Daaaaave: http://github.com/u003f/daaaaave/releases/tag/original
    """
    reader = libsbml.SBMLReader()
    sbml = reader.readSBMLFromFile(filename)
    return sbml


def block_all_imports(sbml):
    """
    Written to mimic neilswainston matlab function blockAllImports
    """
    model = sbml.getModel()

    # strip out format used in recon 2.1
    species = model.getSpecies('M_carbon_e')
    if species:
        species.setBoundaryCondition(True)

    for reaction in model.getListOfReactions():
        nR, nP = 0, 0
        for reactant in reaction.getListOfReactants():
            sID = reactant.getSpecies()
            if not model.getSpecies(sID).getBoundaryCondition():
                nR += 1
        for product in reaction.getListOfProducts():
            sID = product.getSpecies()
            if not model.getSpecies(sID).getBoundaryCondition():
                nP += 1
        kineticLaw = reaction.getKineticLaw()
        if (nR == 1) and (nP == 0):
            kineticLaw.getParameter('LOWER_BOUND').setValue(0)
        if (nR == 0) and (nP == 1):
            kineticLaw.getParameter('UPPER_BOUND').setValue(0)


def change_rxn_bounds(sbml, rxn_name_list, value, bound_type='b'):
    """
    Written to mimic the matlab function changeRxnBounds from http://opencobra.sf.net/
    """
    model = sbml.getModel()
    # convert single entries to lists
    if isinstance(rxn_name_list, str):
        rxn_name_list = [rxn_name_list]
    if isinstance(value, (int, float, long, complex)):
        value = [value] * len(rxn_name_list)
    if isinstance(bound_type, str):
        bound_type = [bound_type] * len(rxn_name_list)
    for index, rID in enumerate(rxn_name_list):
        reaction = get_reaction_by_id(sbml, rID)
        if not reaction:
            print 'reaction %s not found'%rID
        else:
            kineticLaw = reaction.getKineticLaw()
            if bound_type[index] in ['l', 'b']:
                kineticLaw.getParameter('LOWER_BOUND').setValue(value[index])
            if bound_type[index] in ['u', 'b']:
                kineticLaw.getParameter('UPPER_BOUND').setValue(value[index])


def change_objective(sbml, rxn_name_list, objective_coeff=1):
    """
    Written to mimic the matlab function changeObjective from http://opencobra.sf.net/
    """
    model = sbml.getModel()
    for reaction in model.getListOfReactions():
        kineticLaw = reaction.getKineticLaw()
        kineticLaw.getParameter('OBJECTIVE_COEFFICIENT').setValue(0)
    # convert single entries to lists
    if isinstance(rxn_name_list, str):
        rxn_name_list = [rxn_name_list]
    if isinstance(objective_coeff, (int, float, long, complex)):
        objective_coeff = [objective_coeff] * len(rxn_name_list)
    for index, rID in enumerate(rxn_name_list):
        reaction = get_reaction_by_id(sbml, rID)
        if not reaction:
            print 'reaction %s not found'%rID
        else:
            kineticLaw = reaction.getKineticLaw()
            kineticLaw.getParameter('OBJECTIVE_COEFFICIENT').setValue(objective_coeff[index])


def format_for_SBML_ID(txt):
    """
    Written to mimic the matlab function formatForSBMLID from http://opencobra.sf.net/
    """
    txt = 'R_' + txt
    for symbol, replacement in [
            ('-', '_DASH_'),
            ('/', '_FSLASH_'),
            ('\\', '_BSLASH_'),
            ('(', '_LPAREN_'),
            (')', '_RPAREN_'),
            ('[', '_LSQBKT_'),
            (']', '_RSQBKT_'),
            (',', '_COMMA_'),
            ('.', '_PERIOD_'),
            ('\'', '_APOS_'),
            ('&', '&amp'),
            ('<', '&lt'),
            ('>', '&gt'),
            ('"', '&quot')]:
        txt = txt.replace(symbol, replacement)
    return txt


def optimize_cobra_model(sbml):
    """
    Replicate Cobra command optimizeCbModel(model,[],'one').
    Copied from Daaaaave: http://github.com/u003f/daaaaave/releases/tag/original
    """

    bound = INF
    cobra = convert_sbml_to_cobra(sbml, bound)

    N, L, U = cobra['S'], list(cobra['lb']), list(cobra['ub'])
    f, b = list(cobra['c']), list(cobra['b'])
    v_sol, f_opt, conv = easy_lp(f, N, b, L, U, one=False)
    return v_sol, f_opt


def convert_sbml_to_cobra(sbml, bound=INF):
    """
    Get Cobra matrices from SBML model.
    Copied from Daaaaave: http://github.com/u003f/daaaaave/releases/tag/original
    """
    model = sbml.getModel()
    S = sparse.lil_matrix((model.getNumSpecies(), model.getNumReactions()))
    lb, ub, c, b, rev, sIDs = [], [], [], [], [], []
    for species in model.getListOfSpecies():
        sIDs.append(species.getId())
        b.append(0.)
    sIDs = [species.getId() for species in model.getListOfSpecies()]
    for j, reaction in enumerate(model.getListOfReactions()):
        for reactant in reaction.getListOfReactants():
            sID = reactant.getSpecies()
            s = reactant.getStoichiometry()
            if not model.getSpecies(sID).getBoundaryCondition():
                i = sIDs.index(sID)
                S[i, j] = S[i, j] - s
        for product in reaction.getListOfProducts():
            sID = product.getSpecies()
            s = product.getStoichiometry()
            if not model.getSpecies(sID).getBoundaryCondition():
                i = sIDs.index(sID)
                S[i, j] = S[i, j] + s
        kinetic_law = reaction.getKineticLaw()
        rxn_lb = kinetic_law.getParameter('LOWER_BOUND').getValue()
        rxn_ub = kinetic_law.getParameter('UPPER_BOUND').getValue()
        rxn_c = kinetic_law.getParameter('OBJECTIVE_COEFFICIENT').getValue()
        rxn_rev = reaction.getReversible()
        if rxn_lb < -bound:
            rxn_lb = -bound
        if rxn_ub > bound:
            rxn_ub = bound
        if rxn_lb < 0:
            rxn_rev = True
        lb.append(rxn_lb)
        ub.append(rxn_ub)
        c.append(rxn_c)
        rev.append(rxn_rev)
    lb, ub, c, b = numpy.array(lb), numpy.array(ub), numpy.array(c), numpy.array(b)
    rev = numpy.array(rev)
    cobra = {'S': S, 'lb': lb, 'ub': ub, 'c': c, 'b': b, 'rev': rev}
    return cobra


def easy_lp(f, a, b, vlb, vub, one=False):
    '''
    Optimize lp using friends of Gurobi.
    Copied from Daaaaave: http://github.com/u003f/daaaaave/releases/tag/original
    '''

    # create gurobi model
    lp = gurobipy.Model()
    lp.Params.OutputFlag = 0
    lp.Params.FeasibilityTol = 1e-9  # as per Cobra
    lp.Params.OptimalityTol = 1e-9  # as per Cobra
    rows, cols = a.shape
    # add variables to model
    for j in xrange(cols):
        LB = vlb[j]
        if LB == -INF:
            LB = -gurobipy.GRB.INFINITY
        UB = vub[j]
        if UB == INF:
            UB = gurobipy.GRB.INFINITY
        lp.addVar(lb=LB, ub=UB, obj=f[j])
    lp.update()
    lpvars = lp.getVars()
    # iterate over the rows of S adding each row into the model
    S = a.tocsr()
    for i in xrange(rows):
        start = S.indptr[i]
        end = S.indptr[i+1]
        variables = [lpvars[j] for j in S.indices[start:end]]
        coeff = S.data[start:end]
        expr = gurobipy.LinExpr(coeff, variables)
        lp.addConstr(lhs=expr, sense=gurobipy.GRB.EQUAL, rhs=b[i])
    lp.update()
    lp.ModelSense = -1
    lp.optimize()

    v = numpy.empty(len(f))
    v[:] = NAN
    f_opt = NAN
    conv = False
    if lp.Status == gurobipy.GRB.OPTIMAL:
        f_opt = lp.ObjVal
        conv = True
        v = [var.x for var in lp.getVars()]

    # remove model: better memory management?
    del lp

    if conv and one:
        # minimise one norm
        col = sparse.lil_matrix(f)
        a = sparse.vstack([a, f])
        b.append(f_opt)
        f = [0.] * len(f)
        nS, nR = a.shape
        for i in xrange(nR):
            col = sparse.lil_matrix((nS + i, 1))
            a = sparse.hstack([a, col, col])
            row = sparse.lil_matrix((1, nR+2*i+2))
            row[0, i] = 1.
            row[0, nR+2*i] = 1.
            row[0, nR+2*i+1] = -1.
            a = sparse.vstack([a, row])
            vlb.append(0.)
            vlb.append(0.)
            vub.append(INF)
            vub.append(INF)
            f.append(-1.)
            f.append(-1.)
            b.append(0.)
        v_sol = easy_lp(f, a, b, vlb, vub, one=False)[0]
        v = v_sol[:nR]

    if f_opt == -0.0:
        f_opt = 0.0

    return v, f_opt, conv


def get_reaction_by_id(sbml, rID):
    model = sbml.getModel()
    reaction = model.getReaction(rID)
    if not reaction:
        # try cobra replacements
        rID = format_for_SBML_ID(rID)
        reaction = model.getReaction(rID)
    if not reaction:
        # try removing trailing underscore
        if rID[-1] == '_':
            rID = rID[:-1]
        reaction = model.getReaction(rID)
    if not reaction:
        # try adding "_in"
        reaction = model.getReaction(rID + '_in')
    if not reaction:
        # try known alternatives
        rID_map = {
            'R_DM_atp_c': 'R_HKt',  # alternative ATPase
            'R_EX_HC02175_LPAREN_e_RPAREN': 'R_EX_dca_LPAREN_e_RPAREN_'  # alternative C10:0
            }
        if rID in rID_map:
            rID = rID_map[rID]
            reaction = get_reaction_by_id(sbml, rID)
    return reaction


def set_import_bounds(sbml, rxn_name_list, value):
    model = sbml.getModel()
    # convert single entries to lists
    if isinstance(rxn_name_list, str):
        rxn_name_list = [rxn_name_list]
    if isinstance(value, (int, float, long, complex)):
        value = [value] * len(rxn_name_list)
    for index, rID in enumerate(rxn_name_list):
        reaction = get_reaction_by_id(sbml, rID)
        if not reaction:
            print 'reaction %s not found'%rID
        else:
            nR, nP = 0, 0
            for reactant in reaction.getListOfReactants():
                sID = reactant.getSpecies()
                if not model.getSpecies(sID).getBoundaryCondition():
                    nR += 1
            for product in reaction.getListOfProducts():
                sID = product.getSpecies()
                if not model.getSpecies(sID).getBoundaryCondition():
                    nP += 1
            kineticLaw = reaction.getKineticLaw()
            val = abs(value[index])
            if (nR == 0) and (nP == 1):
                kineticLaw.getParameter('UPPER_BOUND').setValue(val)
            elif (nR == 1) and (nP == 0):
                kineticLaw.getParameter('LOWER_BOUND').setValue(-val)
            else:
                print 'reaction %s not import'%rID


def set_infinite_bounds(sbml):
    """Set default bounds to INF, rather than 1000 (say)"""
    model = sbml.getModel()
    for reaction in model.getListOfReactions():
        kineticLaw = reaction.getKineticLaw()
        param = kineticLaw.getParameter('LOWER_BOUND')
        if param.getValue() < -100:
            param.setValue(-INF)
        param = kineticLaw.getParameter('UPPER_BOUND')
        if param.getValue() > 100:
            param.setValue(INF)


if __name__ == '__main__':
    run_all()
    print 'DONE!'

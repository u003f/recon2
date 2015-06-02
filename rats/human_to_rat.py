"""
Create rat model from recon 2.2
Updated: 1 Jun 15 by Kieran Smallbone
"""

from sympy.logic import boolalg
import csv
import numpy as np
import os
import re
import scipy.sparse as sparse
import time

# http://www.gurobi.com/documentation/6.0/quickstart_mac/py_python_interface.html
import gurobipy
# http://frank-fbergmann.blogspot.co.uk/2014/05/libsbml-python-bindings-5101.html
import libsbml


INF = np.inf
NAN = np.nan


def make_rat():
    """
    Create rat model from recon 2.2
    """
    rats_path = os.path.dirname(__file__)
    models_path = os.path.join(rats_path, '..', 'models')
    recon22_path = os.path.join(models_path, 'recon_2.2.xml')
    rambo_path = os.path.join(rats_path, 'rambo.xml')
    sbml = read_sbml(recon22_path)

#     v_sol, f_opt = optimize_cobra_model(sbml, 1000)
#     print 'optimal growth rate:\t%g\n'%(f_opt)

    # map hgnc to rat
    hgnc_map = {}
    hgnc_map_file = os.path.join(rats_path, 'ensembl_hgnc_to_rat.txt')
    with open(hgnc_map_file, 'rU') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            #  ensembl = row[0].strip()
            hgnc = row[1].strip()
            hgnc = hgnc.replace('HGNC:','')
            if len(row) > 2:
                rat = row[2].strip()
            else:
                rat = ''

            if hgnc not in hgnc_map:
                hgnc_map[hgnc] = ''
            if not hgnc_map[hgnc]:
                hgnc_map[hgnc] = rat

    model = sbml.getModel()
    to_remove = []
    for reaction in model.getListOfReactions().clone():
        rID = reaction.getId()
        ga_human = get_notes_field(rID, 'GENE_ASSOCIATION', sbml)
        if ga_human:
            ga_rat = ga_human
            for hgnc in re.findall(r'\b\S+\b', ga_human):
                if (hgnc not in ['and', 'or', 'AND', 'OR']):
                    if hgnc in hgnc_map:
                        rat = hgnc_map[hgnc]
                    else:
                        print 'HGNC:%s\tnot found'%hgnc
                        rat = ''
                    if not rat:
                        rat = 'False'
                    ga_rat = re.sub(r'\b' + hgnc + r'\b', rat, ga_rat)
            ga_rat = to_dnf(ga_rat)
            set_notes_field(rID, 'GENE_ASSOCIATION', ga_rat, sbml)
#             if not ga_rat:
#                 to_remove.append(rID)
#                 print '%s\t[%s]\n%s\n'%(rID, reaction.getName(), ga_human)
#     print '%s of %s reactions removed\n'%(len(to_remove), model.getNumReactions())

    for rID in to_remove:
        model.removeReaction(rID)

    model.setName('rambo')
    model.setId('rambo')
    model.setMetaId('meta_' + model.getId())
    model.setAnnotation(model.getAnnotation())

    write_sbml(sbml, rambo_path)

    v_sol, f_opt = optimize_cobra_model(sbml, 1000)
    print 'optimal growth rate:\t%g\n'%(f_opt)


def to_dnf(association):

    # A and B and (C or D) or E
    association = association.replace(' AND ',' & ').replace(' OR ',' | ').replace(' and ',' & ').replace(' or ',' | ')
    # -> A & B & (C | D) | E
    association = str(boolalg.to_dnf(association))
    # -> Or(And(A, B, C), And(A, B, D), E)
    for and_old in re.findall(r'And\([^)]+\)',association):
        and_new = and_old
        and_new = and_new.replace('And(','(')
        and_new = and_new.replace(', ',' and ')
        association = association.replace(and_old, and_new)
    # -> Or((A and B and C), (A and B and D), E)
    association = association.replace(', ', ' or ')
    if association[:3] == 'Or(':
        association = association[3:-1]
    # .. -> (A and B) or (A and C) or D

    if association == 'False':
        association = ''
    return association


def get_notes_field(eID, name, sbml):
    """
    """
    element = sbml.getModel().getElementBySId(eID)
    try:
        notes = element.getNotesString()
        f = re.search(name + ':([^<]+)', notes)
        f = f.group(1)
        f = f.strip()
    except:
         f = ''
    return f


def set_notes_field(eID, name, value, sbml):

    element = sbml.getModel().getElementBySId(eID)
    notes = element.getNotesString()
    if value:
        value = ' ' + value
    notes = re.sub(name + ':[^<]+', name + ':' + value, notes)
    notes = notes.replace(name + ':' + '</p>', name + ':' + value + '</p>')
    element.setNotes(notes)


def read_sbml(filename):
    """
    Read an SBML file from specified path.
    Copied from Daaaaave: http://github.com/u003f/daaaaave/releases/tag/original
    """
    reader = libsbml.SBMLReader()
    sbml = reader.readSBMLFromFile(filename)
    return sbml


def write_sbml(sbml, filename):
    """
    Write an SBML object to the specified path.
    """
    change_modified_date(sbml)
    writer = libsbml.SBMLWriter()
    writer.writeSBMLToFile(sbml, filename)


def change_modified_date(sbml):
    """
    Change SBML modified date to now.
    """
    history = sbml.getModel().getModelHistory()
    if history:
        history.setModifiedDate(libsbml.Date(w3c_time()))
        # remove all but final modified date
        while history.getListModifiedDates().getSize() > 1:
            history.getListModifiedDates().remove(0)


def w3c_time():
    """
    Return time as a string in W3CDTF UTC format.
    """
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def easy_lp(f, a, b, vlb, vub, one=False):
    '''Optimize lp using friends of Gurobi.'''

    # catch np arrays
    f, b, vlb, vub = list(f), list(b), list(vlb), list(vub)

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

    v = np.empty(len(f))
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

    return v, f_opt, conv


def convert_sbml_to_cobra(sbml, bound=INF):
    """Get Cobra matrices from SBML model."""
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
    lb, ub, c, b = np.array(lb), np.array(ub), np.array(c), np.array(b)
    rev = np.array(rev)
    cobra = {'S': S, 'lb': lb, 'ub': ub, 'c': c, 'b': b, 'rev': rev}
    return cobra


def optimize_cobra_model(sbml, bound=INF):
    """Replicate Cobra command optimizeCbModel(model,[],'one')."""

    cobra = convert_sbml_to_cobra(sbml, bound)

    N, L, U = cobra['S'], list(cobra['lb']), list(cobra['ub'])
    f, b = list(cobra['c']), list(cobra['b'])
    v_sol, f_opt, conv = easy_lp(f, N, b, L, U, one=True)

    return v_sol, f_opt

if __name__ == '__main__':
    make_rat()

    print 'DONE!'

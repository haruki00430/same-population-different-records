import numpy as np

def record_binary(truth, sensitivity, specificity, uniform):
    return uniform < np.where(truth, sensitivity, 1-specificity)

def record_patient_view(latent, se_y, se_c, sp_c, se_d, fpr_d, se_n, u):
    """Latent immutable object is shared, never re-created or modified."""
    return {'latent': latent,
            'outcome': record_binary(latent.outcome, se_y, 1, u['y']),
            'c': record_binary(latent.c, se_c, sp_c, u['c']),
            'e': record_binary(latent.eligibility, se_d, 1-fpr_d, u['e']),
            'a': record_binary(latent.action, se_n, 1, u['a'])}

def bernoulli_pair(a, b, independent=False):
    joint = a*b if independent else min(a, b)
    return [(0,0,1-a-b+joint),(1,0,a-joint),(0,1,b-joint),(1,1,joint)]

def paired_counts(rng, counts, a, b, independent=False):
    """Exact collapsed simulation of shared patients, not a normal approximation."""
    both = a*b if independent else min(a,b)
    # numpy multinomial accepts array-valued n in modern NumPy.
    out = rng.multinomial(counts, np.maximum([both, a-both, b-both, 1-a-b+both], 0))
    return out[...,0]+out[...,1], out[...,0]+out[...,2]

def process_joint(e, qa, qb, dna, dnb, fpa, fpb, nna, nnb, independent=False):
    """3x3 observed category probabilities: outside denominator, D only, N&D.

    Enumerates shared eligibility/action states and two independent uniforms
    for denominator and numerator capture. Clinical action can improve using
    monotone potential actions; only recording differs when qa == qb.
    """
    out = np.zeros((3,3))
    for eligible, pe in [(0,1-e),(1,e)]:
        action = bernoulli_pair(qa,qb) if eligible else [(0,0,1.)]
        for aa, ab, pa in action:
            for da, db, pd in bernoulli_pair(dna if eligible else fpa, dnb if eligible else fpb, independent):
                for na, nb, pn in bernoulli_pair(nna*aa,nnb*ab,independent):
                    ca = 0 if not da else (2 if na else 1)
                    cb = 0 if not db else (2 if nb else 1)
                    out[ca,cb] += pe*pa*pd*pn
    assert abs(out.sum()-1) < 1e-12
    return out

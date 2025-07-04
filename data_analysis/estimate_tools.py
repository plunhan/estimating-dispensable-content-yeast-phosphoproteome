import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import gaussian_kde

def estimate_pi_mixture_model(conditional_scores: list[float], 
							  universal_scores: list[float], 
							  phospho_scores: list[float]) -> float: 
	def neg_log_likelihood(pi):
		pi = pi[0]
		if not (0 <= pi <= 1):
			return np.inf
		kde_cond = gaussian_kde(conditional_scores)
		kde_univ = gaussian_kde(universal_scores)
		mix_pdf = pi * kde_cond(phospho_scores) + (1-pi) * kde_univ(phospho_scores)
		return -np.sum(np.log(mix_pdf + 1e-12))

	result = minimize(neg_log_likelihood, x0=[0.5], bounds=[(0, 1)])
	return result.x[0]
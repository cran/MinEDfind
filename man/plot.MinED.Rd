\name{plot.MinED}
\alias{plot.MinED}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the simulation results for nonparametric two-stage Bayesian adaptive designs
}
\description{
Plot the objects returned by other functions, including (1) operating characteristics of the design, including selection percentage and the number of patients treated at each dose; (2) the estimates of toxicity and response probability for each dose in the admissable set and corresponding 95\% credible interval
}
\usage{
\method{plot}{MinED}(x, name, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
the object returned by other functions
  }
  \item{name}{
the name in the object to be plotted
  }
  \item{...}{
ignored arguments
  }
}

\value{
plot.MinED() returns a figure
}
\references{
Rongji Mu, Guoying Xu, Haitao Pan (2020). A nonparametric two-stage Bayesian adaptive design for minimum effective dose (MinED)-based
dosing-finding trials, (under review)
}
\author{
Chia-Wei Hsu, Fang Wang, Rongji Mu, Haitao Pan, Guoying Xu
}

\examples{
## select the MinED based on the trial data
n = c(3, 6, 0, 0, 0)
y = c(0, 1, 0, 0, 0)
z = c(0, 1, 0, 0, 0)
phi_t = 0.3
phi_e = 0.3
eps_t = 0.1 * phi_t
eps_e = 0.1 * phi_e
select.dose <- select.MinED(n, y, z, phi_t, phi_e, eps_t, eps_e, ct = 0.95)
plot.MinED(select.dose)

## get the operating characteristics for nonparametric two-stage Bayesian adaptive designs
ttox = c(0.05, 0.15, 0.3, 0.45, 0.6)
teff = c(0.05, 0.15, 0.3, 0.45, 0.6)
phi_t = 0.3
phi_e = 0.3
eps_t = 0.1 * phi_t
eps_e = 0.1 * phi_e

oc = get.OC.MinED(ttox = ttox, teff = teff, phi_t = phi_t, phi_e = phi_e,
                  eps_t = eps_t, eps_e = eps_e, cohortsize=3, ncohort1 = 6,
                  ncohort2 = 14, ntrial = 100)

plot.MinED(oc, "Sel\%")
plot.MinED(oc, "#Pts.treated")
plot.MinED(oc, "#Pts.response.to.tox")
plot.MinED(oc, "#Pts.response.to.eff")

}


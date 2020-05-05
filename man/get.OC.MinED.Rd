\name{get.OC.MinED}
\alias{get.OC.MinED}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate operating characteristics for finding the minimum effective dose (MinED)
}
\description{
Obtain the operating characteristics of the nonparametric two-stage Bayesian adaptive design for minimum effective dose (MinED)-based dosing-finding trials
}
\usage{
get.OC.MinED(ttox, teff, phi_t, phi_e, ct = 0.95, eps_t, eps_e, d0 = 1,
             cohortsize = 3, ncohort1, ncohort2, ntrial = 100, extrasafe = TRUE,
             cutoff.eli = 0.95, n.earlystop = 12)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ttox}{
a vector containing the true toxicity rates of the investigational dose levels
  }
  \item{teff}{
a vector containing the true response rates of the investigational dose levels
  }
  \item{phi_t}{
the target DLT rate
  }
  \item{phi_e}{
the target response rate
  }
  \item{ct}{
the cutoff used to eliminate the dose for too toxicity. The default value is ct = 0.95
  }
  \item{eps_t}{
a small value such that (phi_t - eps_t, phi_t + eps_t) is an indifference interval of phi_t. The default value is eps_t = 0.1 * phi_t
  }
  \item{eps_e}{
a small value such that (phi_e - eps_e, phi_e + eps_e) is an indifference interval of phi_e. The default value is eps_e = 0.1 * phi_e
  }
  \item{d0}{
the starting dose level. The default value is d0 = 1
  }
  \item{cohortsize}{
the cohort size
  }
  \item{ncohort1}{
the number of cohort used in stage I
  }
  \item{ncohort2}{
the number of cohort used in stage II
  }
  \item{ntrial}{
the number of simulated trial
  }
  \item{extrasafe}{
extrasafe	set extrasafe = TRUE to impose a more stringent stopping rule
  }
  \item{cutoff.eli}{
the cutoff to eliminate an overly toxic dose for safety. The default value is cutoff.eli = 0.95
  }
  \item{n.earlystop}{
the early stopping parameter. The default value is n.earlystop = 12
  }

}

\value{
\code{get.oc.MinED()} returns the operating characteristics of nonparametric two-stage Bayesian adaptive design as a matrix object,
       including:
       (1) true DLT rate at each dose level,
       (2) true efficacy rate at each dose level,
       (3) selection percentage at each dose level,
       (4) the average number of patients treated at each dose level,
       (5) the average number of patients responded to toxicity at each dose level,
       (6) the average number of patients responded to efficacy at each dose level
}
\references{
Rongji Mu, Guoying Xu, Haitao Pan (2020). A nonparametric two-stage Bayesian adaptive design for minimum effective dose (MinED)-based
dosing-finding trials, (under review)
}
\author{
Chia-Wei Hsu, Fang Wang, Rongji Mu, Haitao Pan, Guoying Xu
}

\examples{
ttox = c(0.05, 0.15, 0.3, 0.45, 0.6)
teff = c(0.05, 0.15, 0.3, 0.45, 0.6)
phi_t = 0.3
phi_e = 0.3
eps_t = 0.1 * phi_t
eps_e = 0.1 * phi_e

oc = get.OC.MinED(ttox = ttox, teff = teff, phi_t = phi_t, phi_e = phi_e,
                  eps_t = eps_t, eps_e = eps_e, cohortsize = 3, ncohort1 = 6,
                  ncohort2 = 14, ntrial = 100)
print(oc)
}


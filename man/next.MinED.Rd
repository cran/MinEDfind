\name{next.MinED}
\alias{next.MinED}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine the dose for the next cohort of new patients for single-agent trials that aim to find a minimum effective dose (MinED)
}
\description{
Determine the dose for the next cohort of new patients for single-agent trials that aim to find a MinED
}
\usage{
next.MinED(n, y, z, d, phi_t, phi_e, eps_t, eps_e, ct = 0.95, N1 = 18)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
a vector of number of patients treated at each dose level
  }
  \item{y}{
a vector of number of patients experiencing the toxicity at each dose level (with the same length as candidate doses)
  }
  \item{z}{
a vector of number of patients showing response at each dose level (with the same length as candidate doses)
  }
  \item{d}{
the starting dose level
  }
  \item{phi_t}{
the target DLT rate
  }
  \item{phi_e}{
the target response rate
  }
  \item{eps_t}{
a small value such that (phi_t - eps_t, phi_t + eps_t) is an indifference interval of phi_t. The default value is eps_t = 0.1 * phi_t
  }
  \item{eps_e}{
a small value such that (phi_e - eps_e, phi_e + eps_e) is an indifference interval of phi_e. The default value is eps_e = 0.1 * phi_e
  }
  \item{ct}{
the cutoff used to eliminate the dose for too toxicity. The default value is ct = 0.95
  }
  \item{N1}{
number of trials in the stage 1. The default value is N1 = 18
  }
}

\value{
next.MinED() returns recommended dose level for the next cohort as a list (\code{$nextdose})
}
\references{
Rongji Mu, Guoying Xu, Haitao Pan (2020). A nonparametric two-stage Bayesian adaptive design for minimum effective dose (MinED)-based
dosing-finding trials, (under review)
}
\author{
Chia-Wei Hsu, Fang Wang, Rongji Mu, Haitao Pan, Guoying Xu
}

\examples{
n = c(3, 6, 0, 0, 0)
y = c(0, 1, 0, 0, 0)
z = c(0, 1, 0, 0, 0)
d = 2
phi_t = 0.3
phi_e = 0.3
eps_t = 0.1 * phi_t
eps_e = 0.1 * phi_e
next.dose <- next.MinED(n = n, y = y, z = z, d = d, phi_t = phi_t,
                        phi_e = phi_e, eps_t = eps_t, eps_e = eps_e)
print(next.dose)
}


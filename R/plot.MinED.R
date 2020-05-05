plot.MinED <- function(x, name, ...){
  Dose_Level <- Lower_Efficacy <- Lower_Toxicity <- Posterior_Efficacy_Est <- NULL
  Posterior_Toxicity_Est <- Sec.. <- Upper_Efficacy <- Upper_Toxicity <- NULL
  X..Pts.response.to.eff <- X..Pts.response.to.tox <- X..Pts.treated <- NULL
  if (class(x)[1] != "MinED"){
    stop("the object putting here is not correct")
  }
  else {
    if (is.matrix(x)){
      df <- data.frame(t(x[, -6]))
      df$Dose_Level <- rownames(df)
      rownames(df) <- c()
      if (tolower(name) == "sel%"){
        p <- ggplot(data = df, aes(x = Dose_Level, y = Sec..)) + geom_bar(stat = "identity") + xlab("Dose Level") + ylab("MinED Selection %")
      }
      else if (tolower(name) == "#pts.treated"){
        p <- ggplot(data = df, aes(x = Dose_Level, y = X..Pts.treated)) + geom_bar(stat = "identity") + xlab("Dose Level") +
             ylab("Number of Patients Treated")
      }
      else if (tolower(name) == "#pts.response.to.tox"){
        p <- ggplot(data = df, aes(x = Dose_Level, y = X..Pts.response.to.tox)) + geom_bar(stat = "identity") + xlab("Dose Level") +
             ylab("Number of Toxicities")
      }
      else if (tolower(name) == "#pts.response.to.eff"){
        p <- ggplot(data = df, aes(x = Dose_Level, y = X..Pts.response.to.eff)) + geom_bar(stat = "identity") + xlab("Dose Level") +
             ylab("Number of Efficacy Responses")
      }
    }
    else if (is.list(x)){
      df <- x[[3]]
      df[, 1] <- factor(df[, 1], levels = 1:length(df[, 1]))
      # make sure for doses which dose not have the trial will be a place holder
      df <- df[rowSums(is.na(df)) == 0, ]
      eff_plot <- ggplot() +
                  geom_errorbar(data = df, aes(x = Dose_Level, ymin = Lower_Efficacy, ymax = Upper_Efficacy), width = 0.2, size = 1, color = "blue") +
                  geom_point(data = df, aes(x = Dose_Level, y = Posterior_Efficacy_Est), size = 4, shape = 21, fill = "white") +
                  ylim(c(0, 1)) + geom_hline(yintercept = x[[2]][, 1], linetype="dashed", color = "red", size=1) +
                  ylab("Posterior Efficacy Est") + xlab("Dose Level") + scale_x_discrete(drop = F)
                                                                      ## make sure to keep those empty doses
      tox_plot <- ggplot() +
                  geom_errorbar(data = df, aes(x = Dose_Level, ymin = Lower_Toxicity, ymax = Upper_Toxicity), width = 0.2, size = 1, color = "blue") +
                  geom_point(data = df, aes(x = Dose_Level, y = Posterior_Toxicity_Est), size = 4, shape = 21, fill = "white") +
                  ylim(c(0, 1)) + geom_hline(yintercept = x[[2]][, 2], linetype = "dashed", color = "red", size = 1) +
                  ylab("Posterior Toxicity Est") + xlab("Dose Level") + scale_x_discrete(drop = F)
                                                                      ## make sure to keep those empty doses
      p <- grid.arrange(eff_plot, tox_plot, nrow = 1)

    }

  }
  p
}

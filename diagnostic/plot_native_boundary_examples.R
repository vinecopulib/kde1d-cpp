arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 4L) {
  stop("usage: script pre-pr.csv bulk.csv expert.csv output.png")
}
library(ggplot2)

estimates <- do.call(rbind, lapply(arguments[1:3], read.csv))
estimates$method <- factor(
  estimates$method,
  levels = c("pre-pr", "bulk", "expert"),
  labels = c("pre-PR log/probit", "new bulk", "new bulk + experts")
)
truth <- unique(estimates[c("support", "scenario", "sample_size", "x", "truth")])

plot <- ggplot(estimates, aes(x, estimate, color = method)) +
  geom_line(linewidth = 0.55) +
  geom_line(data = truth, aes(x = x, y = truth), inherit.aes = FALSE,
            color = "black", linewidth = 0.65, linetype = 2) +
  facet_wrap(~scenario, scales = "free", ncol = 2) +
  scale_color_manual(values = c("#777777", "#377eb8", "#e41a1c")) +
  labs(x = NULL, y = "density", color = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")

ggsave(arguments[[4L]], plot, width = 8, height = 10, dpi = 160)

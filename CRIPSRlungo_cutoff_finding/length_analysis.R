library(MASS)
library(ggplot2)
library(scales)
library(knitr)
library(kableExtra)

threshold <- 10
#Set output directory:
output_directory <- ""


#Import control file:
gang_bao_control_ins <- scan("/gangbao_BCL11A_CONTROL_ins_lengths.txt", what = numeric())

# step 1. subtract off 1 from dataset
y<-gang_bao_control_ins-1

# step 2. fit the NB model
.fit <- glm.nb(y ~ 1)

# step 3. extract model parameters
mu_hat <- exp(coef(.fit)[[1]])
theta_hat <- .fit$theta

# step 4. create plot
x_range <- seq(0, 11)
shifted_nb_df <- data.frame(
  count = x_range,
  expected = dnbinom(x = x_range, size = theta_hat, mu = mu_hat) * length(y)
)

p_model_untrans <- ggplot(data = data.frame(v = y), mapping = aes(x = v)) +
  geom_histogram(binwidth = 1, col = "black", fill = "white") +
  theme_bw() +
  geom_line(
  data = shifted_nb_df,
  aes(x = count, y = expected),
  linewidth = 0.9, col = "firebrick2"
  ) +
  ylab("Count") + xlab("Allele length")

p_model_untrans <- p_model_untrans +
  scale_x_continuous(
    breaks = scales::breaks_width(500), # tick every 1
    labels = function(x) x + 1            # show +1 on the axis
  ) +
  xlab("Allele length")

p_model_trans <- p_model_untrans + scale_y_continuous(transform = pseudo_log_trans())



# step 5. compute the right tail probability at the threshold (keep the minus 1 below; we need to subtract 1 for the shifted NB model)
p_val <- stats::pnbinom(q = threshold - 1, size = theta_hat, mu = mu_hat, lower.tail = FALSE)
print(p_val)

# step 6. output the results
ggsave(filename = paste0(output_directory, "/untrans_hist.png"),
       plot = p_model_untrans, device = "png", scale = 1, height = 3, dpi = 330, width = 4)
ggsave(filename = paste0(output_directory, "/trans_hist.png"),
       plot = p_model_trans, device = "png", scale = 1, height = 3, dpi = 330, width = 4)
write.table(p_val, file = paste0(output_directory, "/probability.txt"))

# step 7. print how many reads are actually greater than threshold and total length dataset
print(sum(y > threshold-1))     # number of values > 10
print(length(y))

# step 8. create a table of possible cutoff allele lengths
allele_len <- 1:11
prob_exact <- dnbinom(allele_len - 1, size = theta_hat, mu = mu_hat)

prob_table <- data.frame(
  allele_length = allele_len,
  probability   = prob_exact
)

# Save it
write.csv(prob_table,
          file = paste0(output_directory, "/nb_exact_length_probs_1_11.csv"),
          row.names = FALSE)
# step 9. output table in visually appearing way

kbl(prob_table,
    col.names = c("Allele length", "Probability"),
    digits = 6, align = "c") |>
  kable_styling(full_width = FALSE)


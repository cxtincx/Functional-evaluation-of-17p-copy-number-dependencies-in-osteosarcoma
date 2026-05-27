##flow 

library(readxl)
library(ggplot2)
data <- read_excel("/Users/catincaapostol/Desktop/Master's/project/catinca flow 15052026/flow MFI 27052026.xls")

library(tidyr)

colnames(data)[1] <- "Sample"

data_clean <- data[!grepl("SD", data$Sample, ignore.case = TRUE), ]

data_long <- pivot_longer(data_clean, 
                          cols = c("T0", "T1"), 
                          names_to = "Timepoint", 
                          values_to = "MFI_Value")
data_long$Sample <- factor(data_long$Sample, levels = unique(data_clean$Sample))

p <- ggplot(data_long, aes(x = Sample, y = MFI_Value, fill = Timepoint)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + 
  scale_fill_manual(values = c("T0" = "forestgreen", "T1" = "grey30")) + 
  theme_minimal() +
  labs(
    x = "Conditions",
    y = "MFI Scale"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 16, color = "black"), 
    legend.position = "top",
    legend.text = element_text (size = 16),
    legend.title = element_text (size = 16),
    axis.title.x = element_text (size = 20), 
    axis.title.y = element_text (size = 20), 
    axis.text.y = element_text (size = 16, color = "black")
  )

print(p)

ggsave(filename = "/Users/catincaapostol/Desktop/Master's/project/flow MFI.png", plot = p, height = 5, width = 10, dpi = 300)

#significance 

library(dplyr)
library(knitr)

sd_row <- data[grepl("sd", data$Sample, ignore.case = TRUE), ]
sd_value <- mean(c(as.numeric(sd_row$t0), as.numeric(sd_row$t1)), na.rm = TRUE)

df_clean <- data[!grepl("sd", data$Sample, ignore.case = TRUE), ]
df_clean$t0 <- as.numeric(df_clean$t0)
df_clean$t1 <- as.numeric(df_clean$t1)

individual_analysis <- df_clean %>%
  mutate(
    Delta = t1 - t0,
    Z_Score = Delta / sd_value,
    Significant = ifelse(abs(Z_Score) >= 1.96, "Yes", "No")
  ) %>%
  dplyr::select(
    `Condition` = Sample,
    `t0 Value` = t0,
    `t1 Value` = t1,
    `Absolute Increase` = Delta,
    `Z-Score` = Z_Score,
    `Significant (p < 0.05)` = Significant
  )

kable(individual_analysis, digits = 4, caption = "Individual Condition Significance Table")
library(dplyr)
library(knitr)

sd_t0 <- 287
sd_t1 <- 255
pooled_sd <- sqrt((sd_t0^2 + sd_t1^2) / 2)

data <- data.frame(
  Sample = c("NRH-OS1 + enAsCas12a", "NRH-OS1 + enAsCas12a negative control", "NRH-OS1 WT", "NRH-OS1 WT negative control"),
  T0 = c(519, 316, 891, 556),
  T1 = c(713, 456, 1158, 696)
)

df_filtered <- data %>%
  mutate(Delta = T1 - T0) %>%
  filter(Sample %in% c("NRH-OS1 + enAsCas12a", "NRH-OS1 WT"))

delta_1 <- df_filtered$Delta[df_filtered$Sample == "NRH-OS1 + enAsCas12a"]
delta_2 <- df_filtered$Delta[df_filtered$Sample == "NRH-OS1 WT"]

diff_increases <- delta_1 - delta_2
se_diff <- sqrt(2) * pooled_sd
z_score <- diff_increases / se_diff
p_val <- 2 * (1 - pnorm(abs(z_score)))

comparison_table <- data.frame(
  `enAsCas12a Increase` = delta_1,
  `WT Increase` = delta_2,
  `Difference` = diff_increases,
  `Z-Score` = z_score,
  `p-value` = p_val,
  `Significant (p < 0.05)` = ifelse(p_val < 0.05, "Yes", "No"),
  check.names = FALSE
)

kable(comparison_table, digits = 4, caption = "Significance Test Between Specific Increases")
write.csv(comparison_table, file = "/Users/catincaapostol/Desktop/Master's/project/mfi stats.csv", row.names = FALSE)




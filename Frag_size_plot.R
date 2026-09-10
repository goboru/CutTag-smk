# workflow/scripts/Frag_size_plot.R

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_pdf  <- args[2]

# -------- CHECK INPUT --------
if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir)
}

# -------- LOAD DATA --------
files <- list.files(
  input_dir,
  pattern = "_fragment_sizes.txt$",
  full.names = TRUE
)

# -------- SAMPLE NAMES --------

sample_names <- basename(files)
sample_names <- sub("_fragment_sizes.txt", "", sample_names)

# -------- MAKE PLOT --------
colors <- rainbow(length(files), alpha = 0.7)

pdf(output_pdf, width = 8, height = 6)

plot(
  NULL,
  xlim = c(0, 1000),
  ylim = c(0, 0.02),
  xlab = "Fragment Size (bp)",
  ylab = "Density",
  main = "ATAC-seq Fragment Size Distribution"
)

for(i in seq_along(files)){

  sizes <- scan(files[i], quiet = TRUE)

  dens <- density(
    sizes,
    from = 0,
    to = 1000,
    bw = 5
  )

  lines(
    dens,
    col = colors[i],
    lwd = 2
  )

  legend(
  	"topright",
  	legend = sample_names,
  	col = colors,
  	lwd = 2,
  	cex = 0.6,
  	ncol = 2,
  	bty = "n"
  )
}

dev.off()




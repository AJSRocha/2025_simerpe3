library(sf)
library(rnaturalearth)
library(rnaturalearthdata)


# desenha o mapa

world <- ne_countries(scale = "medium", returnclass = "sf")
portugal <- subset(world, admin == "Portugal")

# coords <- matrix(c(
#   -9.000, 37.030,  # A
#   -9.282, 36.805,  # B
#   -7.392, 36.805,  # C
#   -7.392, 37.030,  # D
#   -9.000, 37.030  # close polygon
# ), 
# ncol = 2, byrow = TRUE)

coords <- matrix(c(
  -9.000, 37.15,  # A
  -9.282, 36.805,  # B
  -7.392, 36.805,  # C
  -7.392, 37.15,  # D
  -9.000, 37.15  # close polygon
), 
ncol = 2, byrow = TRUE)

# Create an sf polygon
poly <- st_polygon(list(coords)) |> 
  st_sfc(crs = 4326) |>       # WGS84 (lat/lon)
  st_sf(name = "ICES_27.9.a.s.a_boundary", geometry = _)

# Plot the polygon


ggplot() +
  geom_sf(data = portugal, fill = "grey85", color = "black") +
  geom_sf(data = poly, fill = "skyblue", color = "blue", alpha = 0.5) +
  coord_sf(xlim = c(-10, -7), ylim = c(36.5, 37.5)) +
  theme_minimal() +
  labs(
    title = "ICES Subdivision 27.9.a.s.a (Simplified Boundary)")
  )



# cria os dados

ficheiros = list.files('data/GPCPMON/')


nc_path <- "data/GPCPMON/"
nc_name <- ficheiros[1] 
nc_fname <- paste(nc_path, nc_name,  sep="")
dname <- "tmp"  # note: tmp means temperature (not temporary)

# read the netCDF file
tmp_raster <- rast(nc_fname)

# tmp_array <- values(tmp_raster, mat=T)

tmp_df <- as.data.frame(tmp_raster, xy = TRUE)


tmp_sf <- st_as_sf(tmp_df, coords = c("x", "y"), crs = 4326)
inside_sf <- tmp_sf[st_within(tmp_sf, poly, sparse = FALSE), ]
inside_df <- st_drop_geometry(inside_sf)
head(inside_df)

names(inside_df)[grepl("precipitation", names(inside_df))]


chuva = inside_df %>% 
  mutate(ponto = 1:4) %>% 
  select(ponto, names(inside_df)[grepl("precipitation", names(inside_df))]) %>% 
  pivot_longer(-c("ponto")) %>% 
  group_by(name) %>% 
  summarise(precip = sum(value)) %>% 
  mutate(year = rep(1994:2024, each = 12),
         month = rep(1:12, times = 31))

chuva %>% 
ggplot() + 
  geom_line(aes(x = paste(year, month, sep = '_'), 
                y = precip,
                group = 1), col = 'lightblue', size = 1) + 
  theme_bw() + 
  labs(x = "month", y = 'precipitation(m)') + 
  theme(axis.text.x = element_blank())

save(chuva, file = 'data/chuva.rdata')







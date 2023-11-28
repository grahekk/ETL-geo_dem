CREATE OR REPLACE FUNCTION get_geocellids_in_continent(continent_name text)
RETURNS TABLE (geocellid character varying(254))
AS $$
BEGIN
    RETURN QUERY EXECUTE
    'SELECT grid.geocellid
     FROM osm.esa_global_dem_grid AS grid
     JOIN osm.world_continents_boundaries AS na ON ST_Intersects(grid.geometry, na.geometry)
     WHERE na.continent = $1'
    USING continent_name;
END;
$$
LANGUAGE plpgsql;

 -- query below
-- SELECT geocellid
-- FROM osm.esa_global_dem_grid AS grid, osm.world_continents_boundaries AS na 
-- WHERE ST_Intersects(grid.geometry, na.geometry) AND na.continent = 'North America';

 -- and function usage below
-- SELECT * FROM get_geocellids_in_continent('North America');

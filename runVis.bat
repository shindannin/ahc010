call environment.bat
FOR /L %%I IN (10,1,10) DO (
java -jar %JAR% -delay 10 -exec %VS_RELEASE% -seed %%I -vis
)

@REM java -jar SmallPolygonsVis.jar -exec "Release\marathon_main.exe" -seed 2 -vis
@REM java -jar SmallPolygonsVis.jar -manual -debug

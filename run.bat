call environment.bat

FOR /L %%I IN (1,1,10) DO (
	java -jar %JAR% -exec %VS_RELEASE% -seed %%I -novis
)


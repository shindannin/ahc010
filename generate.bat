call environment.bat

cd java
FOR /L %%I IN (1,1,10) DO (
java -jar %JAR% -novis -seed %%I -exec "%VS_RELEASE% %EXAMPLE%\example%%I.txt"
)


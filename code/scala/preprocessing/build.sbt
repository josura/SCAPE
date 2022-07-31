scalaVersion := "2.13.8"
name := "preprocessing"
organization := "unict.it"
version := "1.0"


libraryDependencies += "org.apache.spark" %% "spark-core" % "3.3.0"
libraryDependencies += "org.apache.spark" %% "spark-sql" % "3.3.0"
libraryDependencies +=  "org.apache.hadoop" % "hadoop-client-api" % "3.3.2"
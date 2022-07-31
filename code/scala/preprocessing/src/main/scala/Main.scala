import org.apache.spark.sql.SparkSession
import org.apache.spark.sql.functions.{udf,col,max}
import org.apache.spark.sql.functions.row_number
import org.apache.spark.sql.functions.expr
import org.apache.spark.sql.expressions.Window

object Main extends App {

  val spark:SparkSession = SparkSession
    .builder()
    .appName("Preprocessing scRNA-seq data and annotation datasets")
    .master("local[*]")
    .getOrCreate()
  
  
  import spark.implicits._

  val homoSapAnnotations = spark.read.option("delimiter", "\t").option("header", "true").csv("/run/media/josura/Codigos/dati/tesi/referenceGenomes/human/referenceTranscriptome/Homo_sapiens.GRCh38.107.refseq.tsv")

  val partitionedAnnotations = homoSapAnnotations.groupBy("transcript_stable_id").
    agg(expr("max_by(gene_stable_id, source_identity)").as("gene_stable_id"))//, max("source_identity").as("source_identity"))

  partitionedAnnotations.write.option("header",true).option("delimiter", "\t").csv("/run/media/josura/Codigos/dati/tesi/referenceGenomes/human/referenceTranscriptome/referenceTranscriptomeAnnotation.tsv")

  println("Hello, World!")

  spark.stop()
}
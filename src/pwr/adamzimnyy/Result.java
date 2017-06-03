package pwr.adamzimnyy;

/**
 * Created by adamz on 02.06.2017.
 */
public class Result {

    public static final String FILE_HEADER = "traits\tmethod\tnormal\tsuccess";

    int traits;
    int method;
    int normalization;
    double success;
    int k;

    public Result(int traits, int method, int normalization, double success) {
        this.traits = traits;
        this.method = method;
        this.normalization = normalization;
        this.success = success;
    }

    public Result(int traits, int method, int normalization, double success, int k) {
        this.traits = traits;
        this.method = method;
        this.normalization = normalization;
        this.success = success;
        this.k = k;
    }

    public double getSuccess() {
        return success;
    }

    public void setSuccess(double success) {
        this.success = success;
    }

    public int getTraits() {
        return traits;
    }

    public void setTraits(int traits) {
        this.traits = traits;
    }

    public int getMethod() {
        return method;
    }

    public void setMethod(int method) {
        this.method = method;
    }

    public int getNormalization() {
        return normalization;
    }

    public void setNormalization(int normalization) {
        this.normalization = normalization;
    }

    @Override
    public String toString() {
        return "\ttraits=" + traits +
                "\tmethod=" + method +
                "\tnormalization=" + normalization +
                "\tsuccess=" + success +
                "\tk=" + k;
    }


    public String toFileString() {
        return traits + "\t" + success;
    }
}

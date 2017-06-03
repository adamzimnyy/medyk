package pwr.adamzimnyy;


import java.util.Arrays;

public class DataEntry {
    int id;
    private double[] data = new double[59];
    private int diagnosis;
    public void setData(String[] data) {
        for (int i = 0; i < data.length; i++) {
            this.data[i] = Integer.parseInt(data[i]);
        }
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public double[] getData() {
        return data;
    }

    public void setData(double[] data) {
        this.data = data;
    }

    public int getDiagnosis() {
        return diagnosis;
    }

    public void setDiagnosis(int diagnosis) {
        this.diagnosis = diagnosis;
    }



    @Override
    public String toString() {
        return "DataEntry{" +
                "\n\tdata=" + Arrays.toString(data) + "," +
                 "\n\tdiagnosis=" + diagnosis +
                "\n}";
    }

}




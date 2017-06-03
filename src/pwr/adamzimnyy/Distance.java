package pwr.adamzimnyy;

/**
 * Created by adamz on 03.06.2017.
 */
public class Distance{

    double distance;
    DataEntry point;

    public Distance(double distance, DataEntry point) {
        this.distance = distance;
        this.point = point;
    }

    public double getDistance() {
        return distance;
    }

    public void setDistance(double distance) {
        this.distance = distance;
    }

    public DataEntry getPoint() {
        return point;
    }

    public void setPoint(DataEntry point) {
        this.point = point;
    }
}
